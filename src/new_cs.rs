use std::path::{Path, PathBuf};

use ahash::HashMap;
use bitvec::vec::BitVec;
use color_eyre::Result;

use itertools::Itertools;
use linnet::half_edge::{
    subgraph::{cycle::SignedCycle, Inclusion, OrientedCut, SubGraph, SubGraphOps},
    EdgeData, EdgeId, Flow, Hedge, HedgeGraph, HedgeGraphBuilder, HedgeVec, Involution, NodeIndex,
    Orientation, Parent, TraversalTree,
};
use symbolica::atom::{representation::InlineNum, Atom};

use crate::{
    feyngen::{
        diagram_generator::FeynGen, FeynGenFilter, FeynGenFilters, FeynGenOptions, GenerationType,
        NumeratorAwareGraphGroupingOption, SelfEnergyFilterOptions, SnailFilterOptions,
        TadpolesFilterOptions,
    },
    graph::{BareGraph, DerivedGraphData, Edge, EdgeType, LoopMomentumBasis, Vertex},
    model::Model,
    momentum::{SignOrZero, Signature},
    new_graph::Graph,
    numerator::{NumeratorState, PythonState, UnInit},
    tests_from_pytest::load_generic_model,
    ProcessSettings,
};

struct ProcessDefinition {
    pub initial_pdgs: Vec<i64>, // Do we want a pub type Pdg = i64;?
    pub final_pdgs: Vec<i64>,
}
struct Process<S: NumeratorState = PythonState> {
    pub definition: ProcessDefinition,
    pub collection: ProcessCollection<S>,
}

struct GenerationOptions {
    pub feyngen_options: FeynGenOptions,
    pub numerator_aware_isomorphism_grouping: NumeratorAwareGraphGroupingOption,
    pub filter_self_loop: bool,
    pub graph_prefix: String,
    pub selected_graphs: Option<Vec<String>>,
    pub vetoed_graphs: Option<Vec<String>>,
    pub loop_momentum_bases: Option<HashMap<String, Vec<String>>>,
}

// should produce e+e- -> d d~ g
fn test_process() -> GenerationOptions {
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
    };

    GenerationOptions {
        feyngen_options,
        numerator_aware_isomorphism_grouping: NumeratorAwareGraphGroupingOption::OnlyDetectZeroes,
        filter_self_loop: true,
        graph_prefix: "GL".to_string(),
        selected_graphs: None,
        vetoed_graphs: None,
        loop_momentum_bases: None,
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
            options.numerator_aware_isomorphism_grouping,
            options.filter_self_loop,
            options.graph_prefix,
            options.selected_graphs,
            options.vetoed_graphs,
            options.loop_momentum_bases,
        )?;

        let hedge_collection = bare_collection
            .into_iter()
            .map(|bare_graph| {
                let mut hedge_graph_builder = HedgeGraphBuilder::<Edge, Vertex>::new();

                for vertex in bare_graph.vertices {
                    // do we need to remove the vertices attached to externals?
                    hedge_graph_builder.add_node(vertex.clone());
                }

                for edge in bare_graph.edges {
                    match edge.edge_type {
                        EdgeType::Virtual => {
                            hedge_graph_builder.add_edge(
                                NodeIndex(edge.vertices[0]),
                                NodeIndex(edge.vertices[1]),
                                edge,
                                true,
                            );
                        }
                        EdgeType::Incoming => {
                            hedge_graph_builder.add_external_edge(
                                NodeIndex(edge.vertices[1]),
                                edge,
                                // I am not sure what to put for these two
                                true,
                                Flow::Source,
                            );
                        }
                        EdgeType::Outgoing => {
                            hedge_graph_builder.add_external_edge(
                                NodeIndex(edge.vertices[0]),
                                edge,
                                // I am not sure what to put for these two
                                true,
                                Flow::Source,
                            );
                        }
                    }
                }

                Graph::new(
                    Atom::parse(&bare_graph.overall_factor)
                        .expect("failed to convert overall factor to symbolica atom"),
                    hedge_graph_builder.build(),
                )
            })
            .collect_vec();

        let collection = match options.feyngen_options.generation_type {
            GenerationType::Amplitude => {
                let amplitude = Amplitude {
                    graphs: hedge_collection,
                };

                todo!("convert to PythonState")
            }
            GenerationType::CrossSection => {
                todo!("build representation of cuts")
            }
        };

        Ok(Process {
            definition,
            collection,
        })
    }
}

struct ProcessList<S: NumeratorState = PythonState> {
    pub processes: Vec<Process<S>>,
}

struct ExportSettings {
    root_folder: PathBuf,
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

enum ProcessCollection<S: NumeratorState = PythonState> {
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

struct Amplitude<S: NumeratorState = PythonState> {
    graphs: Vec<Graph<S>>,
}

impl<S: NumeratorState> Amplitude<S> {
    pub fn new() -> Self {
        Self { graphs: vec![] }
    }

    pub fn add_graph(&mut self, graph: Graph<S>) -> Result<()> {
        self.graphs.push(graph);
        /// TODO: validate that the graph is compatible
        Ok(())
    }
}
struct CrossSection<S: NumeratorState = PythonState> {
    supergraphs: Vec<CrossSectionGraph<S>>,
}

pub struct CrossSectionGraph<S: NumeratorState = PythonState> {
    graph: Graph<S>,
    cuts: Vec<CrossSectionCut>,
}

pub struct CrossSectionCut {
    cut: OrientedCut,
    left: BitVec,
    right: BitVec,
}

#[test]
fn test_new_structure() {
    let model = load_generic_model("sm");
    let generation_options = test_process();

    let mut process_list = ProcessList::new();
    process_list.generate(generation_options, &model).unwrap();
}
