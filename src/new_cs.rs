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
    numerator::{NumeratorState, PythonState, UnInit},
    tests_from_pytest::load_generic_model,
    ProcessSettings,
};

pub struct Graph<S: NumeratorState = PythonState> {
    multiplicity: Atom,
    underlying: HedgeGraph<Edge, Vertex>,
    loop_momentum_basis: HedgeLMB,
    derived_data: DerivedGraphData<S>,
}

pub trait FeynmanGraph {
    fn new_lmb(&self) -> HedgeLMB;
}

pub struct HedgeLMB {
    tree: TraversalTree,
    lmb_basis: Vec<Hedge>,
}

impl FeynmanGraph for HedgeGraph<Edge, Vertex> {
    fn new_lmb(&self) -> HedgeLMB {
        // root node should contain a dangling (external edge), that will be the dependent external
        let root = self
            .iter_nodes()
            .find(|(a, _)| {
                self.iter_egdes(*a)
                    .any(|(e, _)| matches!(e, EdgeId::Unpaired { .. }))
            })
            .unwrap_or(self.iter_nodes().next().unwrap())
            .0;

        let tree = TraversalTree::dfs(self, &self.full_filter(), root, None);

        let leaves = tree.leaf_edges();
        let mut external_leaves = self.empty_filter();
        let mut externals = vec![];

        for i in leaves.included_iter() {
            if self.involution.is_identity(i) {
                external_leaves.set(i.0, true);
                externals.push(i);
            }
        }

        let mut ext_signatures: HedgeVec<Signature> = self.new_derived_edge_data_empty();

        let empty_signature = Signature::from_iter(externals.iter().map(|e| SignOrZero::Zero));
        for (i, &h) in externals.iter().enumerate() {
            let mut signature = empty_signature.clone();
            signature.0[i] = SignOrZero::Plus;
            ext_signatures[h] = Some(signature);
        }

        let mut current_leaf_nodes = tree.leaf_nodes(&self);
        while let Some(leaf_node) = current_leaf_nodes.pop() {
            let hairs = &self.hairs_from_id(leaf_node).hairs;
            let mut root_pointer = None;
            let mut signature = empty_signature.clone();

            for h in hairs.included_iter() {
                match tree.parent(h) {
                    Parent::Root => {}
                    Parent::Hedge { hedge_to_root, .. } => {
                        if *hedge_to_root != h {
                            signature.sum(&ext_signatures[h].as_ref().unwrap());
                        } else {
                            if self
                                .involved_node_hairs(h)
                                .unwrap()
                                .hairs
                                .included_iter()
                                .all(|a| a != h && ext_signatures.is_set(a))
                            {
                                current_leaf_nodes.push(self.involved_node_id(h).unwrap())
                            }
                            root_pointer = Some(h);
                        }
                    }
                    Parent::Unset => {
                        panic!("Unset parent")
                    }
                }
            }

            if let Some(root_pointer) = root_pointer {
                ext_signatures[root_pointer] = Some(signature);
            }
        }

        let tree_complement = tree.tree.complement(self);

        let mut cycle_basis = Vec::new();
        let mut lmb_basis = Vec::new();

        for (e, d) in self.iter_egdes(&tree_complement) {
            if let EdgeId::Paired { source, sink } = e {
                lmb_basis.push(source);
                cycle_basis.push(
                    SignedCycle::from_cycle(tree.cycle(source).unwrap(), source, self).unwrap(),
                );
            }
        }

        let signatures = self
            .involution
            .map_data_ref(|a| (), &|e| ())
            .map_edge_data(|e, d| {
                let e = match e {
                    EdgeId::Paired { source, sink } => {
                        let mut internal = vec![];
                        for (i, c) in cycle_basis.iter().enumerate() {
                            if c.filter.includes(&source) {
                                internal.push(SignOrZero::Plus);
                            } else if c.filter.includes(&sink) {
                                internal.push(SignOrZero::Minus);
                            } else {
                                internal.push(SignOrZero::Zero);
                            }
                        }

                        let internal_signature = Signature::from_iter(internal);

                        // return EdgeData::new(Signature::from_iter(iter), orientation)
                    }
                    EdgeId::Unpaired { hedge, flow } => {}
                    _ => {}
                };
                d
            });

        HedgeLMB { tree, lmb_basis }
    }
}

impl Graph<UnInit> {
    pub fn new(multiplicity: Atom, underlying: HedgeGraph<Edge, Vertex>) -> Self {
        Self {
            multiplicity,
            loop_momentum_basis: underlying.new_lmb(),
            underlying,
            derived_data: DerivedGraphData::new_empty(),
        }
    }
}

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

                hedge_graph_builder.build()
            })
            .collect_vec();

        let collection = match options.feyngen_options.generation_type {
            GenerationType::Amplitude => {
                ProcessCollection::Amplitudes(todo!("turn baregraph into amplitude graph"))
            }
            GenerationType::CrossSection => {
                ProcessCollection::CrossSections(todo!("turn baregraph into crosssectiongraph"))
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
