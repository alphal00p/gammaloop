use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    fmt,
    hash::{Hash, Hasher},
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    },
};

use feynkit_graph::{
    DiagramEdge, DiagramError, DiagramVertex, ExternalState, FeynmanDiagram, ParticleReference,
};
use feynkit_model::{Model, ModelError, Particle, VertexRule};
use linnet::half_edge::involution::Flow;
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike};
use linnet::half_edge::{HedgeGraph, NodeIndex};
use rayon::{ThreadPoolBuildError, ThreadPoolBuilder, prelude::*};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    graph::{GenerationSettings, Graph, HalfEdge},
    parser::ParseSettings,
    symbol,
};
use thiserror::Error;

use crate::{
    FilterScope, GenerationControl, GenerationFilter, GenerationFilterKind, GenerationOptions,
    GenerationProgress, GenerationType, GroupingError, NumeratorGrouping, ParticleSelector,
    Process, ProcessError, SelfEnergyFilterOptions, SewnFilterOptions, SnailFilterOptions,
    TadpoleFilterOptions, grouping,
};

#[derive(Debug, Error)]
pub enum GenerationError {
    #[error(transparent)]
    Process(#[from] ProcessError),
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error(transparent)]
    Diagram(#[from] DiagramError),
    #[error("failed to build generation thread pool: {0}")]
    ThreadPool(#[from] ThreadPoolBuildError),
    #[error("graph automorphism factor cannot be represented by u64: {0}")]
    SymmetryFactor(String),
    #[error("no graph was generated for the requested process")]
    NoGraphs,
    #[error("could not assign an interaction to node {node} with signature {signature:?}")]
    MissingInteraction { node: usize, signature: Vec<String> },
    #[error("option {option} for filter {filter} in the {scope} scope is not implemented")]
    UnsupportedFilterOption {
        scope: FilterScope,
        filter: GenerationFilterKind,
        option: &'static str,
    },
    #[error("filter {filter} is not valid in the {scope} scope for {generation_type:?} generation")]
    InvalidFilterScope {
        scope: FilterScope,
        filter: GenerationFilterKind,
        generation_type: GenerationType,
    },
    #[error("duplicate filter {filter} in the {scope} scope")]
    DuplicateFilter {
        scope: FilterScope,
        filter: GenerationFilterKind,
    },
    #[error(
        "filter {filter} in the {scope} scope has an invalid range: {minimum} is greater than {maximum}"
    )]
    InvalidFilterRange {
        scope: FilterScope,
        filter: GenerationFilterKind,
        minimum: usize,
        maximum: usize,
    },
    #[error("generation thread count must be greater than zero, received {0}")]
    InvalidThreadCount(usize),
    #[error("numerator grouping requires at least one numerical sample")]
    InvalidNumericalSampleCount,
    #[error("arithmetic overflow while computing {0}")]
    ArithmeticOverflow(&'static str),
    #[error(transparent)]
    Grouping(#[from] GroupingError),
    #[error(
        "fermion flow branches at vertex {vertex} with degree {degree}; branching fermion chains are not supported"
    )]
    UnsupportedFermionBranching { vertex: usize, degree: usize },
    #[error("spanning-tree edge between nodes {node} and {parent} is missing")]
    MissingSpanningTreeEdge { node: usize, parent: usize },
    #[error("external node {node} has no incident edge")]
    MissingExternalEdge { node: usize },
    #[error("external leg {index} has no resolved particle")]
    MissingExternalParticle { index: usize },
    #[error("external fermion flow must have exactly two legs, found {legs:?}")]
    InvalidExternalFermionLegCount { legs: Vec<usize> },
    #[error(
        "external fermion flow between legs {legs:?} must connect a fermion and an antifermion, found particles {particles:?}"
    )]
    InvalidExternalFermionPair {
        legs: [usize; 2],
        particles: [i64; 2],
    },
    #[error("sewn external fermion leg {leg} has degree {degree}, expected two")]
    InvalidExternalFermionLoopDegree { leg: usize, degree: usize },
    #[error("failed to reconstruct the normalized fermion-flow graph: {0}")]
    GraphConstruction(String),
    #[error("failed to parse the {owner} numerator '{expression}': {message}")]
    NumeratorParse {
        owner: String,
        expression: String,
        message: String,
    },
    #[error("the {owner} numerator references UFO leg {leg}, but it has {legs} legs")]
    InvalidNumeratorLeg {
        owner: String,
        leg: i64,
        legs: usize,
    },
    #[error("the {owner} numerator contains an invalid UFO index '{index}'")]
    InvalidNumeratorIndex { owner: String, index: String },
    #[error("the {owner} numerator references momentum but has no incident edge")]
    MissingNumeratorMomentum { owner: String },
    #[error(
        "could not assign particle '{particle}' at UFO leg {leg} of interaction '{interaction}' on vertex {vertex}"
    )]
    MissingInteractionLeg {
        vertex: usize,
        interaction: String,
        leg: usize,
        particle: String,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GroupMember {
    /// Stable generated-order index assigned before zero-numerator removal.
    pub source_diagram: usize,
    /// Index into [`GenerationResult::diagrams`] after zero-numerator removal.
    pub diagram: usize,
    /// The numerator of this member divided by the numerator of the master.
    pub ratio: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DiagramGroup {
    /// Post-zero-removal index of the numerator reference diagram.
    pub master: usize,
    /// Deterministically ordered members, including the master with ratio one.
    pub members: Vec<GroupMember>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenerationReport {
    pub topology_count: usize,
    pub interaction_assignment_count: usize,
    /// Number of diagrams retained after exact zero-numerator removal.
    pub retained_count: usize,
    pub zero_numerator_count: usize,
    pub completed: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenerationResult {
    /// Every retained diagram in generated order; grouping does not collapse this vector.
    pub diagrams: Vec<FeynmanDiagram>,
    pub groups: Vec<DiagramGroup>,
    pub report: GenerationReport,
}

#[derive(Debug, Clone, Copy)]
struct NumeratorHalfEdge {
    edge: usize,
    flow: Flow,
}

#[derive(Debug, Clone, Copy)]
enum NumeratorOwner {
    Vertex(usize),
    Edge(usize),
}

impl fmt::Display for NumeratorOwner {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Vertex(vertex) => write!(formatter, "vertex {vertex}"),
            Self::Edge(edge) => write!(formatter, "edge {edge}"),
        }
    }
}

struct NumeratorInstantiation<'a> {
    owner: NumeratorOwner,
    legs: &'a [NumeratorHalfEdge],
}

impl NumeratorInstantiation<'_> {
    // UFO `idx(shift, leg)` and `dummy(index)` values are local to one
    // propagator or interaction. Localize them before exposing the annotation
    // so equal endpoint indices contract and unrelated dummy indices cannot.
    fn instantiate(&self, expression: &str) -> Result<String, GenerationError> {
        let owner = self.owner.to_string();
        let mut atom = Atom::parse(
            expression,
            "feynkit_generator_numerator",
            ParseSettings::default(),
        )
        .map_err(|error| GenerationError::NumeratorParse {
            owner: owner.clone(),
            expression: expression.to_owned(),
            message: error.to_string(),
        })?;
        let mut error = None;
        let mut changed = false;
        atom = atom.replace_map(|term, _, out| {
            if error.is_some() {
                return;
            }
            match self.localize_momentum(term) {
                Ok(Some(replacement)) => {
                    changed = true;
                    **out = replacement;
                }
                Ok(None) => {}
                Err(localization_error) => error = Some(localization_error),
            }
        });
        if let Some(error) = error {
            return Err(error);
        }
        atom = atom.replace_map(|term, _, out| {
            if error.is_some() {
                return;
            }
            match self.localize_indices(term) {
                Ok(Some(replacement)) => {
                    changed = true;
                    **out = replacement;
                }
                Ok(None) => {}
                Err(localization_error) => error = Some(localization_error),
            }
        });
        if let Some(error) = error {
            Err(error)
        } else if changed {
            // The default `Display` printer hides every namespace. Keep the
            // public FeynKit and UFO namespaces in the serialized annotation.
            Ok(atom.to_plain_string())
        } else {
            // Keep already-global model expressions byte-for-byte stable.
            Ok(expression.to_owned())
        }
    }

    fn localize_momentum(&self, term: AtomView<'_>) -> Result<Option<Atom>, GenerationError> {
        match term {
            AtomView::Fun(function) if Self::is_model_symbol(function.get_symbol(), "P") => {
                let arguments = function.iter().collect::<Vec<_>>();
                let momentum = match arguments.as_slice() {
                    [] => self.momentum(self.default_leg()?, None)?,
                    [index] => self.momentum(
                        self.default_leg()?,
                        Some(self.localize_index_argument(*index)?),
                    )?,
                    [index, momentum_leg] => {
                        let leg = self.leg(self.local_leg_number(*momentum_leg)?)?;
                        let momentum =
                            self.momentum(leg, Some(self.localize_index_argument(*index)?))?;
                        // The UFO all-incoming convention reverses momentum on
                        // the source half-edge of an interaction.
                        if leg.flow == Flow::Source {
                            -momentum
                        } else {
                            momentum
                        }
                    }
                    _ => return Err(self.invalid_index(term)),
                };
                Ok(Some(momentum))
            }
            AtomView::Fun(function)
                if Self::is_model_symbol(function.get_symbol(), "PSlash") =>
            {
                let arguments = function.iter().collect::<Vec<_>>();
                if !(2..=3).contains(&arguments.len()) {
                    return Err(self.invalid_index(term));
                }
                let leg = if let Some(momentum_leg) = arguments.get(2) {
                    self.leg(self.local_leg_number(*momentum_leg)?)?
                } else {
                    self.default_leg()?
                };
                let mut localized = FunctionBuilder::new(function.get_symbol());
                for index in &arguments[..2] {
                    localized = localized.add_arg(self.localize_index_argument(*index)?);
                }
                localized = localized.add_arg(self.momentum(leg, None)?);
                Ok(Some(localized.finish()))
            }
            AtomView::Fun(function)
                if Self::is_model_symbol(function.get_symbol(), "Slash")
                    && function.get_nargs() == 1
                    && function.iter().next().is_some_and(|argument| {
                        matches!(argument, AtomView::Var(variable) if Self::is_model_symbol(variable.get_symbol(), "P"))
                    }) =>
            {
                Ok(Some(
                    FunctionBuilder::new(function.get_symbol())
                        .add_arg(self.momentum(self.default_leg()?, None)?)
                        .finish(),
                ))
            }
            _ => Ok(None),
        }
    }

    fn localize_indices(&self, term: AtomView<'_>) -> Result<Option<Atom>, GenerationError> {
        let AtomView::Fun(function) = term else {
            return Ok(None);
        };
        if Self::is_model_symbol(function.get_symbol(), "idx") {
            return self.localize_index_function(term).map(Some);
        }
        if Self::is_model_symbol(function.get_symbol(), "dummy") {
            return self.localize_dummy_function(term).map(Some);
        }
        let symbol = function.get_symbol();
        let name = symbol.get_stripped_name();
        if !Self::is_model_symbol(symbol, name)
            || !matches!(
                name,
                "Identity"
                    | "IdentityL"
                    | "Gamma"
                    | "Gamma5"
                    | "ProjM"
                    | "ProjP"
                    | "Sigma"
                    | "C"
                    | "Metric"
                    | "Epsilon"
                    | "T"
                    | "f"
                    | "d"
                    | "EpsilonBar"
                    | "T6"
                    | "K6"
                    | "K6Bar"
            )
        {
            return Ok(None);
        }

        let mut localized = FunctionBuilder::new(function.get_symbol());
        for argument in function.iter() {
            localized = localized.add_arg(self.localize_index_argument(argument)?);
        }
        Ok(Some(localized.finish()))
    }

    fn localize_index_argument(&self, argument: AtomView<'_>) -> Result<Atom, GenerationError> {
        if let Ok(index) = i64::try_from(argument) {
            return if index > 0 {
                self.index(self.leg(index)?, 1)
            } else if index < 0 {
                self.dummy(
                    usize::try_from(index.unsigned_abs()).map_err(|_| {
                        GenerationError::ArithmeticOverflow("a symbolic dummy index")
                    })?,
                )
            } else {
                Err(self.invalid_index(argument))
            };
        }
        if let AtomView::Fun(function) = argument {
            if Self::is_model_symbol(function.get_symbol(), "idx") {
                return self.localize_index_function(argument);
            }
            if Self::is_model_symbol(function.get_symbol(), "dummy") {
                return self.localize_dummy_function(argument);
            }
        }
        Ok(argument.to_owned())
    }

    fn localize_index_function(&self, index: AtomView<'_>) -> Result<Atom, GenerationError> {
        let AtomView::Fun(function) = index else {
            return Err(self.invalid_index(index));
        };
        let arguments = function.iter().collect::<Vec<_>>();
        let [shift, leg] = arguments.as_slice() else {
            return Err(self.invalid_index(index));
        };
        let shift = i64::try_from(*shift).map_err(|_| self.invalid_index(index))?;
        if shift <= 0 {
            return Err(self.invalid_index(index));
        }
        self.index(
            self.leg(i64::try_from(*leg).map_err(|_| self.invalid_index(index))?)?,
            shift,
        )
    }

    fn localize_dummy_function(&self, dummy: AtomView<'_>) -> Result<Atom, GenerationError> {
        let AtomView::Fun(function) = dummy else {
            return Err(self.invalid_index(dummy));
        };
        let arguments = function.iter().collect::<Vec<_>>();
        let [index] = arguments.as_slice() else {
            return Err(self.invalid_index(dummy));
        };
        let index = i64::try_from(*index).map_err(|_| self.invalid_index(dummy))?;
        if index <= 0 {
            return Err(self.invalid_index(dummy));
        }
        self.dummy(
            usize::try_from(index)
                .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic dummy index"))?,
        )
    }

    fn local_leg_number(&self, argument: AtomView<'_>) -> Result<i64, GenerationError> {
        if let Ok(leg) = i64::try_from(argument) {
            return Ok(leg);
        }
        let AtomView::Fun(function) = argument else {
            return Err(self.invalid_index(argument));
        };
        if !Self::is_model_symbol(function.get_symbol(), "idx") {
            return Err(self.invalid_index(argument));
        }
        let arguments = function.iter().collect::<Vec<_>>();
        let [_, leg] = arguments.as_slice() else {
            return Err(self.invalid_index(argument));
        };
        i64::try_from(*leg).map_err(|_| self.invalid_index(argument))
    }

    fn leg(&self, one_based: i64) -> Result<NumeratorHalfEdge, GenerationError> {
        if one_based > 0 {
            let index = usize::try_from(one_based)
                .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic UFO leg"))?;
            if let Some(leg) = self.legs.get(index - 1) {
                return Ok(*leg);
            }
        }
        Err(GenerationError::InvalidNumeratorLeg {
            owner: self.owner.to_string(),
            leg: one_based,
            legs: self.legs.len(),
        })
    }

    fn default_leg(&self) -> Result<NumeratorHalfEdge, GenerationError> {
        self.legs
            .first()
            .copied()
            .ok_or_else(|| GenerationError::MissingNumeratorMomentum {
                owner: self.owner.to_string(),
            })
    }

    fn index(&self, leg: NumeratorHalfEdge, shift: i64) -> Result<Atom, GenerationError> {
        let edge = i64::try_from(leg.edge)
            .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic edge identifier"))?;
        let symbol = match leg.flow {
            Flow::Source => symbol!("FeynKit::SourceIndex"),
            Flow::Sink => symbol!("FeynKit::SinkIndex"),
        };
        Ok(symbol.call((edge, shift)))
    }

    fn dummy(&self, local: usize) -> Result<Atom, GenerationError> {
        let (symbol, owner) = match self.owner {
            NumeratorOwner::Vertex(vertex) => (symbol!("FeynKit::VertexDummy"), vertex),
            NumeratorOwner::Edge(edge) => (symbol!("FeynKit::EdgeDummy"), edge),
        };
        let owner = i64::try_from(owner)
            .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic numerator owner"))?;
        let local = i64::try_from(local)
            .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic dummy index"))?;
        Ok(symbol.call((owner, local)))
    }

    fn momentum(
        &self,
        leg: NumeratorHalfEdge,
        index: Option<Atom>,
    ) -> Result<Atom, GenerationError> {
        let edge = i64::try_from(leg.edge)
            .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic edge identifier"))?;
        let mut momentum = FunctionBuilder::new(symbol!("FeynKit::Momentum")).add_arg(edge);
        if let Some(index) = index {
            momentum = momentum.add_arg(index);
        }
        Ok(momentum.finish())
    }

    fn invalid_index(&self, index: AtomView<'_>) -> GenerationError {
        GenerationError::InvalidNumeratorIndex {
            owner: self.owner.to_string(),
            index: index.to_string(),
        }
    }

    fn is_model_symbol(symbol: symbolica::atom::Symbol, name: &str) -> bool {
        symbol.get_stripped_name() == name
            && (symbol.get_namespace() == "feynkit_generator_numerator"
                || symbol.get_namespace() == "UFO"
                || symbol.get_namespace().starts_with("UFO::"))
    }
}

#[derive(Clone)]
pub struct Generator {
    model: Arc<Model>,
}

impl Generator {
    pub fn new(model: impl Into<Arc<Model>>) -> Self {
        Self {
            model: model.into(),
        }
    }

    pub fn model(&self) -> &Model {
        &self.model
    }

    pub fn generate(
        &self,
        process: &Process,
        options: &GenerationOptions,
    ) -> Result<GenerationResult, GenerationError> {
        process.validate()?;
        self.validate_options(process, options)?;
        let resolved = ResolvedProcess::new(&self.model, process)?;
        let signatures = InteractionSignatures::new(&self.model, options)?;
        let external_edges = resolved.external_edges(&self.model)?;
        let generation_signatures: Vec<_> = signatures.keys().cloned().collect();

        let progress_count = Arc::new(AtomicUsize::new(0));
        let callback = options.progress.clone();
        let cancellation = options.cancellation.clone();
        let progress_counter = progress_count.clone();
        let progress_cancellation = cancellation.clone();
        let progress_cancellation_check = options.cancellation_check.clone();
        let mut settings = GenerationSettings::new()
            .max_loops(*process.loop_count().end())
            .allow_self_loops(options.allow_self_loops)
            .allow_zero_flow_edges(options.allow_zero_flow_edges)
            .progress_fn(Box::new(move |_| {
                let generated_graphs = progress_counter.fetch_add(1, Ordering::Relaxed) + 1;
                let control = callback
                    .as_ref()
                    .map_or(GenerationControl::Continue, |callback| {
                        callback(GenerationProgress { generated_graphs })
                    });
                if control == GenerationControl::Cancel {
                    progress_cancellation.cancel();
                }
                if progress_cancellation_check
                    .as_ref()
                    .is_some_and(|check| check())
                {
                    progress_cancellation.cancel();
                }
                !progress_cancellation.is_cancelled()
            }));
        if let Some(maximum) = options.max_vertices {
            // Symbolica counts the degree-one external nodes as vertices. The
            // public FeynKit option counts interaction vertices instead.
            settings = settings.max_vertices(maximum.checked_add(external_edges.len()).ok_or(
                GenerationError::ArithmeticOverflow("the Symbolica vertex limit"),
            )?);
        }
        if let Some(maximum) = options.max_bridges() {
            settings = settings.max_bridges(maximum);
        }
        let abort_cancellation = cancellation.clone();
        let abort_cancellation_check = options.cancellation_check.clone();
        settings = settings.abort_check(Box::new(move || {
            abort_cancellation.is_cancelled()
                || abort_cancellation_check
                    .as_ref()
                    .is_some_and(|check| check())
        }));

        let generated = Graph::generate(&external_edges, &generation_signatures, settings);
        let (raw_graphs, mut completed) = match generated {
            Ok(graphs) => (graphs, true),
            Err(graphs) => (graphs, false),
        };
        if options.cancellation_requested() {
            completed = false;
        }
        if raw_graphs.is_empty() && completed {
            return Err(GenerationError::NoGraphs);
        }

        let mut raw_graphs: Vec<_> = raw_graphs.into_iter().collect();
        // Sort for reproducible generated IDs, then discard contributions below
        // the requested loop range before more expensive processing.
        raw_graphs.sort_by(|left, right| left.0.cmp(&right.0));
        raw_graphs.retain(|(graph, _)| process.loop_count().contains(&graph.num_loops()));
        let mut filtered_graphs = Vec::with_capacity(raw_graphs.len());
        for entry in raw_graphs {
            if options.cancellation_requested() {
                completed = false;
                break;
            }
            if self.passes_topology_filters(&entry.0, options)? {
                filtered_graphs.push(entry);
            }
        }
        let raw_graphs = filtered_graphs;
        let topology_count = raw_graphs.len();

        let mut pool = ThreadPoolBuilder::new();
        // Respect an explicit thread count; otherwise Rayon uses its default
        // logical-CPU count.
        if let Some(threads) = options.threads {
            pool = pool.num_threads(threads);
        }
        let pool = pool.build()?;
        let colored = pool.install(|| {
            raw_graphs
                .par_iter()
                .filter_map(|(graph, symmetry)| {
                    if options.cancellation_requested() {
                        return None;
                    }
                    let symmetry = match u64::try_from(symmetry.clone()) {
                        Ok(symmetry) => symmetry,
                        Err(error) => {
                            return Some(Err(GenerationError::SymmetryFactor(error.to_string())));
                        }
                    };
                    Some(
                        self.assign_interactions(graph, &signatures)
                            .map(|graphs| (graphs, symmetry)),
                    )
                })
                .collect::<Result<Vec<_>, GenerationError>>()
        })?;
        if options.cancellation_requested() {
            completed = false;
        }
        let mut colored: Vec<_> = colored
            .into_iter()
            .flat_map(|(graphs, symmetry)| {
                graphs
                    .into_iter()
                    .map(move |(graph, multiplicity)| ColoredTopology {
                        graph,
                        symmetry,
                        multiplicity,
                        cut_partitions: Vec::new(),
                    })
            })
            .collect();
        colored.sort_by(|left, right| left.graph.cmp(&right.graph));
        let interaction_assignment_count = colored.len();

        let mut diagram_inputs = Vec::new();
        for mut topology in colored {
            if options.cancellation_requested() {
                completed = false;
                break;
            }
            if process.generation_type() == GenerationType::CrossSection {
                topology.cut_partitions =
                    resolved.cut_partitions(&topology.graph, &self.model, options)?;
                if topology.cut_partitions.is_empty()
                    || !self.retain_cut_amplitude_partitions(&mut topology, options)?
                    || !resolved.passes_sewn_filter(&topology.graph, options)
                {
                    continue;
                }
            }
            let fermion_loop_count = self.closed_fermion_loop_count(&topology.graph)?;
            if !self.passes_interaction_filters(&topology.graph, fermion_loop_count, options)? {
                continue;
            }
            let normalized = resolved.normalize_fermion_flows(&topology.graph, &self.model)?;
            topology.graph = normalized.graph;
            diagram_inputs.push((topology, fermion_loop_count, normalized.signs));
        }
        let width = diagram_inputs.len().saturating_sub(1).to_string().len();
        let mut diagrams = Vec::with_capacity(diagram_inputs.len());
        for (index, (topology, fermion_loop_count, signs)) in diagram_inputs.into_iter().enumerate()
        {
            if options.cancellation_requested() {
                completed = false;
                break;
            }
            diagrams.push(self.to_diagram(
                format!("{}{index:0width$}", options.graph_prefix),
                topology,
                fermion_loop_count,
                signs,
            )?);
        }
        let grouped = if options.cancellation_requested() {
            completed = false;
            grouping::group_diagrams(diagrams, &self.model, &NumeratorGrouping::None)?
        } else {
            grouping::group_diagrams(diagrams, &self.model, &options.numerator_grouping)?
        };
        if options.cancellation_requested() {
            completed = false;
        }

        Ok(GenerationResult {
            report: GenerationReport {
                topology_count,
                interaction_assignment_count,
                retained_count: grouped.diagrams.len(),
                zero_numerator_count: grouped.zero_numerator_count,
                completed,
            },
            diagrams: grouped.diagrams,
            groups: grouped.groups,
        })
    }

    fn validate_options(
        &self,
        process: &Process,
        options: &GenerationOptions,
    ) -> Result<(), GenerationError> {
        if options.threads == Some(0) {
            return Err(GenerationError::InvalidThreadCount(0));
        }
        if matches!(
            &options.numerator_grouping,
            NumeratorGrouping::Identical(options)
                | NumeratorGrouping::UpToSign(options)
                | NumeratorGrouping::UpToScalar(options)
                if options.number_of_numerical_samples == 0
        ) {
            return Err(GenerationError::InvalidNumericalSampleCount);
        }
        self.validate_filters(
            FilterScope::Graph,
            process.generation_type(),
            &options.graph_filters,
        )?;
        self.validate_filters(
            FilterScope::CutAmplitude,
            process.generation_type(),
            &options.cut_amplitude_filters,
        )
    }

    fn validate_filters(
        &self,
        scope: FilterScope,
        generation_type: GenerationType,
        filters: &[GenerationFilter],
    ) -> Result<(), GenerationError> {
        let mut seen = BTreeSet::new();
        for filter in filters {
            let kind = filter.kind();
            if !seen.insert(kind) {
                return Err(GenerationError::DuplicateFilter {
                    scope,
                    filter: kind,
                });
            }

            match scope {
                FilterScope::Graph => match filter {
                    GenerationFilter::SelfEnergy(options) => {
                        if options.only_scaleless {
                            return Err(GenerationError::UnsupportedFilterOption {
                                scope,
                                filter: kind,
                                option: "only_scaleless",
                            });
                        }
                    }
                    GenerationFilter::Tadpoles(options) => {
                        if options.only_scaleless {
                            return Err(GenerationError::UnsupportedFilterOption {
                                scope,
                                filter: kind,
                                option: "only_scaleless",
                            });
                        }
                    }
                    GenerationFilter::ZeroSnails(options) => {
                        if options.only_scaleless {
                            return Err(GenerationError::UnsupportedFilterOption {
                                scope,
                                filter: kind,
                                option: "only_scaleless",
                            });
                        }
                    }
                    GenerationFilter::ParticleVeto(_)
                    | GenerationFilter::VertexAllow(_)
                    | GenerationFilter::VertexVeto(_)
                    | GenerationFilter::MaxNumberOfBridges(_)
                    | GenerationFilter::CouplingOrders(_)
                    | GenerationFilter::LoopCountRange(_)
                    | GenerationFilter::FermionLoopCountRange(_)
                    | GenerationFilter::FactorizedLoopTopologiesCountRange(_) => {}
                    GenerationFilter::BlobRange(_) | GenerationFilter::SpectatorRange(_)
                        if generation_type == GenerationType::CrossSection => {}
                    GenerationFilter::BlobRange(_) | GenerationFilter::SpectatorRange(_) => {
                        return Err(GenerationError::InvalidFilterScope {
                            scope,
                            filter: kind,
                            generation_type,
                        });
                    }
                    GenerationFilter::PerturbativeOrders(_) | GenerationFilter::Sewn(_)
                        if generation_type == GenerationType::CrossSection => {}
                    GenerationFilter::PerturbativeOrders(_) | GenerationFilter::Sewn(_) => {
                        return Err(GenerationError::InvalidFilterScope {
                            scope,
                            filter: kind,
                            generation_type,
                        });
                    }
                },
                FilterScope::CutAmplitude => {
                    if generation_type != GenerationType::CrossSection {
                        return Err(GenerationError::InvalidFilterScope {
                            scope,
                            filter: kind,
                            generation_type,
                        });
                    }
                    match filter {
                        GenerationFilter::CouplingOrders(_)
                        | GenerationFilter::LoopCountRange(_) => {}
                        GenerationFilter::SelfEnergy(_)
                        | GenerationFilter::Tadpoles(_)
                        | GenerationFilter::ZeroSnails(_)
                        | GenerationFilter::Sewn(_)
                        | GenerationFilter::ParticleVeto(_)
                        | GenerationFilter::VertexAllow(_)
                        | GenerationFilter::VertexVeto(_)
                        | GenerationFilter::MaxNumberOfBridges(_)
                        | GenerationFilter::BlobRange(_)
                        | GenerationFilter::SpectatorRange(_)
                        | GenerationFilter::PerturbativeOrders(_)
                        | GenerationFilter::FermionLoopCountRange(_)
                        | GenerationFilter::FactorizedLoopTopologiesCountRange(_) => {
                            return Err(GenerationError::InvalidFilterScope {
                                scope,
                                filter: kind,
                                generation_type,
                            });
                        }
                    }
                }
            }

            self.validate_filter_range(scope, filter)?;
            match filter {
                GenerationFilter::ParticleVeto(pdgs) => {
                    for pdg in pdgs {
                        self.model.particle_by_pdg(*pdg)?;
                    }
                }
                GenerationFilter::VertexAllow(names) | GenerationFilter::VertexVeto(names) => {
                    for name in names {
                        self.model.vertex_rule(name)?;
                    }
                }
                GenerationFilter::CouplingOrders(required) => {
                    for name in required.keys() {
                        self.model.order(name)?;
                    }
                }
                GenerationFilter::PerturbativeOrders(required) => {
                    for name in required.keys() {
                        self.model.order(name)?;
                    }
                }
                GenerationFilter::SelfEnergy(_)
                | GenerationFilter::Tadpoles(_)
                | GenerationFilter::ZeroSnails(_)
                | GenerationFilter::Sewn(_)
                | GenerationFilter::MaxNumberOfBridges(_)
                | GenerationFilter::LoopCountRange(_)
                | GenerationFilter::BlobRange(_)
                | GenerationFilter::SpectatorRange(_)
                | GenerationFilter::FermionLoopCountRange(_)
                | GenerationFilter::FactorizedLoopTopologiesCountRange(_) => {}
            }
        }
        Ok(())
    }

    fn validate_filter_range(
        &self,
        scope: FilterScope,
        filter: &GenerationFilter,
    ) -> Result<(), GenerationError> {
        let invalid = match filter {
            GenerationFilter::CouplingOrders(required) => {
                required
                    .values()
                    .find_map(|(minimum, maximum)| match maximum {
                        Some(maximum) if minimum > maximum => Some((*minimum, *maximum)),
                        _ => None,
                    })
            }
            GenerationFilter::LoopCountRange((minimum, maximum))
            | GenerationFilter::FermionLoopCountRange((minimum, maximum))
            | GenerationFilter::FactorizedLoopTopologiesCountRange((minimum, maximum)) => {
                (minimum > maximum).then_some((*minimum, *maximum))
            }
            GenerationFilter::BlobRange(range) | GenerationFilter::SpectatorRange(range) => {
                (range.start() > range.end()).then_some((*range.start(), *range.end()))
            }
            GenerationFilter::SelfEnergy(_)
            | GenerationFilter::Tadpoles(_)
            | GenerationFilter::ZeroSnails(_)
            | GenerationFilter::Sewn(_)
            | GenerationFilter::ParticleVeto(_)
            | GenerationFilter::VertexAllow(_)
            | GenerationFilter::VertexVeto(_)
            | GenerationFilter::MaxNumberOfBridges(_)
            | GenerationFilter::PerturbativeOrders(_) => None,
        };
        if let Some((minimum, maximum)) = invalid {
            return Err(GenerationError::InvalidFilterRange {
                scope,
                filter: filter.kind(),
                minimum,
                maximum,
            });
        }
        Ok(())
    }

    fn assign_interactions(
        &self,
        graph: &Graph<NodeColor, EdgeColor>,
        signatures: &InteractionSignatures,
    ) -> Result<Vec<InteractionAssignment>, GenerationError> {
        let mut assignments: Vec<Vec<ColoredNode>> = vec![Vec::new()];
        for (node_index, node) in graph.nodes().iter().enumerate() {
            let choices = if let Some(external) = &node.data.external {
                vec![ColoredNode::External(external.clone())]
            } else {
                let mut signature = Vec::new();
                for edge_id in &node.edges {
                    let edge = &graph.edges()[*edge_id];
                    let direction = edge.directed.then_some(edge.vertices.0 == node_index);
                    signature.push(EdgeColor {
                        pdg: edge.data.pdg,
                        direction,
                    });
                    // A self-loop contributes both half-edges to the
                    // interaction signature.
                    if edge.vertices.0 == edge.vertices.1 {
                        signature.push(EdgeColor {
                            pdg: edge.data.pdg,
                            direction: direction.map(|value| !value),
                        });
                    }
                }
                signature.sort();
                signatures
                    .rules(&signature)
                    .ok_or_else(|| GenerationError::MissingInteraction {
                        node: node_index,
                        signature: signature
                            .iter()
                            .map(|edge| format!("{:?}:{}", edge.direction, edge.pdg))
                            .collect(),
                    })?
                    .iter()
                    .cloned()
                    .map(ColoredNode::Interaction)
                    .collect()
            };

            assignments = assignments
                .into_iter()
                .flat_map(|assignment| {
                    choices.iter().cloned().map(move |choice| {
                        let mut next = assignment.clone();
                        next.push(choice);
                        next
                    })
                })
                .collect();
        }

        let mut grouped = BTreeMap::new();
        for assignment in assignments {
            let mut colored = Graph::new();
            for node in assignment {
                colored.add_node(node);
            }
            for edge in graph.edges() {
                colored
                    .add_edge(edge.vertices.0, edge.vertices.1, edge.directed, edge.data)
                    .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
            }
            *grouped.entry(colored.canonize().graph).or_insert(0) += 1;
        }
        Ok(grouped.into_iter().collect())
    }

    fn passes_topology_filters(
        &self,
        graph: &Graph<NodeColor, EdgeColor>,
        options: &GenerationOptions,
    ) -> Result<bool, GenerationError> {
        let self_energy = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::SelfEnergy(options) => Some(options),
                _ => None,
            });
        let tadpoles = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::Tadpoles(options) => Some(options),
                _ => None,
            });
        let zero_snails = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::ZeroSnails(options) => Some(options),
                _ => None,
            });
        let factorized_loop_range = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::FactorizedLoopTopologiesCountRange(range) => Some(range),
                _ => None,
            });
        if self.veto_special_topologies(
            graph,
            self_energy,
            tadpoles,
            zero_snails,
            factorized_loop_range,
        )? {
            return Ok(false);
        }

        for filter in &options.graph_filters {
            match filter {
                GenerationFilter::ParticleVeto(vetoed)
                    if graph
                        .edges()
                        .iter()
                        .any(|edge| vetoed.contains(&edge.data.pdg)) =>
                {
                    return Ok(false);
                }
                GenerationFilter::LoopCountRange((minimum, maximum))
                    if !(*minimum..=*maximum).contains(&graph.num_loops()) =>
                {
                    return Ok(false);
                }
                GenerationFilter::ParticleVeto(_)
                | GenerationFilter::VertexAllow(_)
                | GenerationFilter::VertexVeto(_)
                | GenerationFilter::MaxNumberOfBridges(_)
                | GenerationFilter::CouplingOrders(_)
                | GenerationFilter::LoopCountRange(_)
                | GenerationFilter::FermionLoopCountRange(_)
                | GenerationFilter::SelfEnergy(_)
                | GenerationFilter::Tadpoles(_)
                | GenerationFilter::ZeroSnails(_)
                | GenerationFilter::Sewn(_)
                | GenerationFilter::BlobRange(_)
                | GenerationFilter::SpectatorRange(_)
                | GenerationFilter::FactorizedLoopTopologiesCountRange(_)
                | GenerationFilter::PerturbativeOrders(_) => {}
            }
        }
        Ok(true)
    }

    fn veto_special_topologies(
        &self,
        graph: &Graph<NodeColor, EdgeColor>,
        veto_self_energy: Option<&SelfEnergyFilterOptions>,
        veto_tadpole: Option<&TadpoleFilterOptions>,
        veto_snails: Option<&SnailFilterOptions>,
        factorized_loop_topologies_count_range: Option<&(usize, usize)>,
    ) -> Result<bool, GenerationError> {
        if veto_self_energy.is_none()
            && veto_tadpole.is_none()
            && veto_snails.is_none()
            && factorized_loop_topologies_count_range.is_none()
        {
            return Ok(false);
        }

        let external_nodes: Vec<_> = graph
            .nodes()
            .iter()
            .enumerate()
            .filter_map(|(index, node)| node.data.external.as_ref().map(|_| index))
            .collect();
        if graph.nodes().is_empty() {
            return Ok(factorized_loop_topologies_count_range
                .is_some_and(|range| !(range.0..=range.1).contains(&0)));
        }
        let roots: Vec<_> = if external_nodes.is_empty() {
            vec![0]
        } else {
            external_nodes
        };
        let mut external_particles = BTreeMap::new();
        for (node_index, node) in graph.nodes().iter().enumerate() {
            if let Some(external) = &node.data.external {
                let edge_id = node
                    .edges
                    .first()
                    .ok_or(GenerationError::MissingExternalEdge { node: node_index })?;
                let edge = &graph.edges()[*edge_id];
                external_particles
                    .insert(external.index, self.model.particle_by_pdg(edge.data.pdg)?);
            }
        }

        // Test vetoing from all external spanning-tree root positions to make
        // sure the result does not depend on spanning-tree directions.
        // TODO rewrite and improve the vetoing logic of special topologies.
        for root in roots {
            if !self.veto_special_topologies_with_spanning_tree_root(
                graph,
                veto_self_energy,
                veto_tadpole,
                veto_snails,
                factorized_loop_topologies_count_range,
                &external_particles,
                root,
            )? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    #[allow(clippy::too_many_arguments)]
    fn veto_special_topologies_with_spanning_tree_root(
        &self,
        graph: &Graph<NodeColor, EdgeColor>,
        veto_self_energy: Option<&SelfEnergyFilterOptions>,
        veto_tadpole: Option<&TadpoleFilterOptions>,
        veto_snails: Option<&SnailFilterOptions>,
        factorized_loop_topologies_count_range: Option<&(usize, usize)>,
        external_particles: &BTreeMap<usize, &Particle>,
        spanning_tree_root: usize,
    ) -> Result<bool, GenerationError> {
        let max_external = external_particles.len();
        let spanning_tree_node_external_index = graph.nodes()[spanning_tree_root]
            .data
            .external
            .as_ref()
            .map(|external| external.index);
        let mut spanning_tree = graph.get_spanning_tree(spanning_tree_root);
        spanning_tree.chain_decomposition();

        let mut node_children = vec![vec![]; spanning_tree.nodes.len()];
        for (node_index, node) in spanning_tree.nodes.iter().enumerate() {
            node_children[node.parent].push(node_index);
        }

        let mut external_momenta_routing = vec![vec![]; spanning_tree.nodes.len()];
        for (node_index, node) in graph.nodes().iter().enumerate() {
            let Some(external) = &node.data.external else {
                continue;
            };
            if node_index == spanning_tree_root {
                continue;
            }
            external_momenta_routing[node_index].push(external.index);
            let mut next_node = spanning_tree.nodes[node_index].parent;
            while next_node != spanning_tree_root {
                external_momenta_routing[next_node].push(external.index);
                next_node = spanning_tree.nodes[next_node].parent;
            }
        }
        if let Some(root_external_index) = spanning_tree_node_external_index {
            for route in &mut external_momenta_routing {
                if route.len() == max_external - 1 {
                    *route = vec![root_external_index];
                }
            }
        }

        // See https://arxiv.org/pdf/1209.0700 for information on the logic of
        // this algorithm.

        // Tuple format: (external leg id, back-edge start node position,
        // back-edge position in list, chain id).
        let mut self_energy_attachments: HashSet<(usize, usize, usize, usize)> = HashSet::new();
        // Tuple format: (back-edge start node position, back-edge position in
        // list, chain id).
        let mut vacuum_attachments: HashSet<(usize, usize, usize)> = HashSet::new();
        // Tuple format: (back-edge start node position, back-edge position in
        // list, chain id).
        let mut self_loops: HashSet<(usize, usize, usize)> = HashSet::new();
        let mut factorized_loop_count = 0;

        let mut visited_nodes = vec![false; spanning_tree.nodes.len()];
        for &node_index in &spanning_tree.order {
            let node = &spanning_tree.nodes[node_index];
            for (back_edge_index, back_edge) in node.back_edges.iter().enumerate() {
                let chain_id = node_index;
                if back_edge.target == node_index {
                    factorized_loop_count += 1;
                    self_loops.insert((node_index, back_edge_index, chain_id));
                    continue;
                }
                let mut self_energy_external_leg_id = None;
                let mut current_chain_node = back_edge.target;
                let mut is_valid_chain = true;
                loop {
                    if current_chain_node == node_index {
                        factorized_loop_count += 1;
                        break;
                    }
                    let momenta = &external_momenta_routing[current_chain_node];
                    if momenta.len() == 1 {
                        if self_energy_external_leg_id.is_some_and(|leg_id| leg_id != momenta[0]) {
                            is_valid_chain = false;
                            break;
                        }
                        self_energy_external_leg_id = Some(momenta[0]);
                        for child in &node_children[current_chain_node] {
                            if !external_momenta_routing[*child].is_empty()
                                && external_momenta_routing[*child][0] != momenta[0]
                            {
                                is_valid_chain = false;
                                break;
                            }
                        }
                        if !is_valid_chain {
                            break;
                        }
                    } else if !momenta.is_empty() {
                        is_valid_chain = false;
                        break;
                    }
                    if spanning_tree.nodes[current_chain_node].chain_id.is_none()
                        || visited_nodes[current_chain_node]
                    {
                        is_valid_chain = false;
                        break;
                    }
                    visited_nodes[current_chain_node] = true;
                    current_chain_node = spanning_tree.nodes[current_chain_node].parent;
                }

                if is_valid_chain {
                    if let Some(leg_id) = self_energy_external_leg_id {
                        // Make sure the self-energy attachment point does not
                        // receive any other external momenta. For 1 -> 1
                        // processes, also verify that it is not the whole graph.
                        if external_momenta_routing[node_index].as_slice() == [leg_id]
                            && (max_external > 2
                                || spanning_tree
                                    .nodes
                                    .iter()
                                    .filter(|node| !node.external && node.chain_id.is_none())
                                    .count()
                                    > 1)
                        {
                            self_energy_attachments.insert((
                                leg_id,
                                node_index,
                                back_edge_index,
                                chain_id,
                            ));
                        }
                    } else {
                        vacuum_attachments.insert((node_index, back_edge_index, chain_id));
                    }
                }
            }
        }

        if factorized_loop_topologies_count_range.is_some_and(|(minimum, maximum)| {
            factorized_loop_count < *minimum || factorized_loop_count > *maximum
        }) {
            return Ok(true);
        }

        let mut tree_bridge_node_indices = HashSet::new();
        for (node_index, node) in spanning_tree.nodes.iter().enumerate() {
            if node.chain_id.is_none()
                && !node.external
                && !external_momenta_routing[node_index].is_empty()
                && !spanning_tree.nodes[node.parent]
                    .back_edges
                    .iter()
                    .any(|end| node_index == end.target)
            {
                tree_bridge_node_indices.insert(node_index);
            }
        }

        // Confirm self-energies by checking that the back-edge start node is a
        // bridge.
        for (leg_id, back_edge_start_node_index, _, _) in &self_energy_attachments {
            if tree_bridge_node_indices.contains(back_edge_start_node_index)
                && let Some(options) = veto_self_energy
            {
                let particle = external_particles
                    .get(leg_id)
                    .ok_or(GenerationError::MissingExternalParticle { index: *leg_id })?;
                if (!particle.is_massless() && options.veto_massive)
                    || (particle.is_massless() && options.veto_massless)
                {
                    return Ok(true);
                }
            }
        }

        // For vacuum attachments, differentiate snails (the start node is a
        // tree bridge) from tadpoles (the start node is not a tree bridge).
        for (back_edge_start_node_index, _, _) in vacuum_attachments.iter().chain(self_loops.iter())
        {
            let attachment_particle_is_massive = if max_external > 0 {
                let mut first_tree_attachment_node_index = *back_edge_start_node_index;
                while external_momenta_routing[first_tree_attachment_node_index].is_empty()
                    && spanning_tree.nodes[first_tree_attachment_node_index]
                        .chain_id
                        .is_none()
                {
                    first_tree_attachment_node_index =
                        spanning_tree.nodes[first_tree_attachment_node_index].parent;
                }
                let parent = spanning_tree.nodes[first_tree_attachment_node_index].parent;
                let edge = graph
                    .edges()
                    .iter()
                    .find(|edge| {
                        edge.vertices == (first_tree_attachment_node_index, parent)
                            || edge.vertices == (parent, first_tree_attachment_node_index)
                    })
                    .ok_or(GenerationError::MissingSpanningTreeEdge {
                        node: first_tree_attachment_node_index,
                        parent,
                    })?;
                !self.model.particle_by_pdg(edge.data.pdg)?.is_massless()
            } else {
                // Always consider the attachment particle massive for vacuum
                // graphs: without external legs there is no supported way to
                // differentiate massive and massless attachments.
                true
            };

            if !tree_bridge_node_indices.contains(back_edge_start_node_index)
                && spanning_tree.nodes[*back_edge_start_node_index]
                    .chain_id
                    .is_none()
            {
                if let Some(options) = veto_tadpole
                    && ((attachment_particle_is_massive && options.veto_attached_to_massive)
                        || (!attachment_particle_is_massive && options.veto_attached_to_massless))
                {
                    return Ok(true);
                }
            } else if let Some(options) = veto_snails
                && ((attachment_particle_is_massive && options.veto_attached_to_massive)
                    || (!attachment_particle_is_massive && options.veto_attached_to_massless))
            {
                return Ok(true);
            }
        }

        Ok(false)
    }

    fn passes_interaction_filters(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
        fermion_loop_count: usize,
        options: &GenerationOptions,
    ) -> Result<bool, GenerationError> {
        for filter in &options.graph_filters {
            match filter {
                GenerationFilter::CouplingOrders(required) => {
                    let actual = self.coupling_orders(graph)?;
                    if !required.iter().all(|(name, (minimum, maximum))| {
                        let value = actual.get(name).copied().unwrap_or(0);
                        value >= *minimum && maximum.is_none_or(|maximum| value <= maximum)
                    }) {
                        return Ok(false);
                    }
                }
                GenerationFilter::FermionLoopCountRange((minimum, maximum)) => {
                    if !(*minimum..=*maximum).contains(&fermion_loop_count) {
                        return Ok(false);
                    }
                }
                GenerationFilter::ParticleVeto(_)
                | GenerationFilter::VertexAllow(_)
                | GenerationFilter::VertexVeto(_)
                | GenerationFilter::MaxNumberOfBridges(_)
                | GenerationFilter::LoopCountRange(_)
                | GenerationFilter::SelfEnergy(_)
                | GenerationFilter::Tadpoles(_)
                | GenerationFilter::ZeroSnails(_)
                | GenerationFilter::Sewn(_)
                | GenerationFilter::BlobRange(_)
                | GenerationFilter::SpectatorRange(_)
                | GenerationFilter::FactorizedLoopTopologiesCountRange(_)
                | GenerationFilter::PerturbativeOrders(_) => {}
            }
        }
        Ok(true)
    }

    fn retain_cut_amplitude_partitions(
        &self,
        topology: &mut ColoredTopology,
        options: &GenerationOptions,
    ) -> Result<bool, GenerationError> {
        let mut retained = Vec::with_capacity(topology.cut_partitions.len());
        for partition in topology.cut_partitions.drain(..) {
            let mut passes = true;
            for filter in &options.cut_amplitude_filters {
                passes &= match filter {
                    GenerationFilter::CouplingOrders(required) => {
                        [&partition.left, &partition.right].iter().all(|side| {
                            required.iter().all(|(name, (minimum, maximum))| {
                                let value = side.coupling_orders.get(name).copied().unwrap_or(0);
                                value >= *minimum && maximum.is_none_or(|maximum| value <= maximum)
                            })
                        })
                    }
                    // The amplitude loop filter constrains the sum of the loop
                    // counts on both sides of the cut, matching GammaLoop's
                    // cross-section generation semantics.
                    GenerationFilter::LoopCountRange((minimum, maximum)) => (*minimum..=*maximum)
                        .contains(&(partition.left.loop_count + partition.right.loop_count)),
                    _ => {
                        return Err(GenerationError::InvalidFilterScope {
                            scope: FilterScope::CutAmplitude,
                            filter: filter.kind(),
                            generation_type: GenerationType::CrossSection,
                        });
                    }
                };
            }
            if passes {
                retained.push(partition);
            }
        }
        topology.cut_partitions = retained;
        Ok(!topology.cut_partitions.is_empty())
    }

    fn coupling_orders(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
    ) -> Result<BTreeMap<String, usize>, GenerationError> {
        let mut result = BTreeMap::new();
        for node in graph.nodes() {
            let ColoredNode::Interaction(rule_name) = &node.data else {
                continue;
            };
            let rule = self.model.vertex_rule(rule_name)?;
            for (name, order) in rule.coupling_orders(&self.model)? {
                *result.entry(name).or_insert(0) += order;
            }
        }
        Ok(result)
    }

    // Count closed, non-branching fermion components. Branching components are
    // rejected explicitly below.
    fn closed_fermion_loop_count(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
    ) -> Result<usize, GenerationError> {
        let mut adjacency = vec![Vec::new(); graph.nodes().len()];
        let mut degrees = vec![0; graph.nodes().len()];
        // Build adjacency using only fermion edges. A self-loop contributes
        // degree two.
        for edge in graph.edges() {
            if !self.model.particle_by_pdg(edge.data.pdg)?.is_fermion() {
                continue;
            }
            let (left, right) = edge.vertices;
            if left == right {
                degrees[left] += 2;
                adjacency[left].push(right);
            } else {
                degrees[left] += 1;
                degrees[right] += 1;
                adjacency[left].push(right);
                adjacency[right].push(left);
            }
        }

        for (vertex, degree) in degrees.iter().copied().enumerate() {
            if degree > 2 {
                return Err(GenerationError::UnsupportedFermionBranching { vertex, degree });
            }
        }

        let mut visited = vec![false; graph.nodes().len()];
        let mut closed_loops = 0;
        for start in 0..graph.nodes().len() {
            if degrees[start] == 0 || visited[start] {
                continue;
            }
            let mut stack = vec![start];
            let mut closed = true;
            while let Some(vertex) = stack.pop() {
                if visited[vertex] {
                    continue;
                }
                visited[vertex] = true;
                closed &= degrees[vertex] == 2;
                stack.extend(
                    adjacency[vertex]
                        .iter()
                        .copied()
                        .filter(|neighbor| !visited[*neighbor]),
                );
            }
            closed_loops += usize::from(closed);
        }
        Ok(closed_loops)
    }

    fn to_diagram(
        &self,
        name: String,
        topology: ColoredTopology,
        fermion_loop_count: usize,
        fermion_signs: FermionSigns,
    ) -> Result<FeynmanDiagram, GenerationError> {
        // Symbolica's automorphism factor is a denominator in the overall
        // multiplier.
        let mut overall_factor = if topology.multiplicity == 1 {
            format!("1/AutG({})", topology.symmetry)
        } else {
            format!(
                "CouplingsMultiplicity({})/AutG({})",
                topology.multiplicity, topology.symmetry
            )
        };
        if fermion_loop_count % 2 == 1 {
            overall_factor = format!("InternalFermionLoopSign(-1)*{overall_factor}");
        }
        if fermion_signs.include_external_ordering {
            let sign = if fermion_signs.external_ordering_negative {
                -1
            } else {
                1
            };
            overall_factor = format!("ExternalFermionOrderingSign({sign})*{overall_factor}");
        }
        if fermion_signs.include_antifermion_spin_sum {
            let sign = if fermion_signs.antifermion_spin_sum_negative {
                -1
            } else {
                1
            };
            overall_factor = format!("AntiFermionSpinSumSign({sign})*{overall_factor}");
        }
        let mut builder = FeynmanDiagram::builder(name)
            .symmetry_factor(topology.symmetry)
            .overall_factor(overall_factor);
        let mut numerator_factors = Vec::new();
        for (index, node) in topology.graph.nodes().iter().enumerate() {
            let vertex = match &node.data {
                ColoredNode::External(external) => DiagramVertex::external(
                    format!("ext{}", external.index),
                    external.index,
                    external.state,
                ),
                ColoredNode::Interaction(rule_name) => {
                    let rule = self.model.vertex_rule(rule_name)?;
                    let raw_numerator = self.interaction_numerator(rule)?;
                    let legs = self.interaction_half_edges(&topology.graph, index, rule)?;
                    let numerator = NumeratorInstantiation {
                        owner: NumeratorOwner::Vertex(index),
                        legs: &legs,
                    }
                    .instantiate(&raw_numerator)?;
                    if numerator != "1" {
                        numerator_factors.push(format!("({numerator})"));
                    }
                    let mut vertex =
                        DiagramVertex::interaction(format!("v{index}"), rule_name.clone());
                    vertex.numerator = Some(numerator);
                    vertex
                }
            };
            builder.add_vertex(vertex);
        }
        for (edge_index, edge) in topology.graph.edges().iter().enumerate() {
            let particle = self.model.particle_by_pdg(edge.data.pdg)?;
            let internal = matches!(
                &topology.graph.nodes()[edge.vertices.0].data,
                ColoredNode::Interaction(_)
            ) && matches!(
                &topology.graph.nodes()[edge.vertices.1].data,
                ColoredNode::Interaction(_)
            );
            let edge_numerator = if internal {
                let legs = [
                    NumeratorHalfEdge {
                        edge: edge_index,
                        flow: Flow::Source,
                    },
                    NumeratorHalfEdge {
                        edge: edge_index,
                        flow: Flow::Sink,
                    },
                ];
                NumeratorInstantiation {
                    owner: NumeratorOwner::Edge(edge_index),
                    legs: &legs,
                }
                .instantiate(self.propagator_numerator(particle)?)?
            } else {
                "1".to_owned()
            };
            if edge_numerator != "1" {
                numerator_factors.push(format!("({edge_numerator})"));
            }
            let mut diagram_edge = DiagramEdge::new(
                ParticleReference::new(&particle.name, particle.pdg_code),
                edge.directed,
            );
            diagram_edge.numerator = Some(edge_numerator);
            builder.add_edge(
                feynkit_graph::VertexId(edge.vertices.0),
                feynkit_graph::VertexId(edge.vertices.1),
                diagram_edge,
            )?;
        }
        let numerator = if numerator_factors.is_empty() {
            "1".to_owned()
        } else {
            numerator_factors.join("*")
        };
        let diagram = builder.numerator(numerator).build()?;
        diagram.validate(&self.model)?;
        Ok(diagram)
    }

    fn interaction_half_edges(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
        vertex: usize,
        rule: &VertexRule,
    ) -> Result<Vec<NumeratorHalfEdge>, GenerationError> {
        // Reconstruct the UFO particle order discarded by topology
        // canonicalization. Direction distinguishes particle and antiparticle
        // slots, while a self-loop contributes both of its half-edges.
        let mut available = Vec::new();
        for edge_id in &graph.nodes()[vertex].edges {
            let edge = &graph.edges()[*edge_id];
            if edge.vertices.0 == vertex {
                available.push((
                    NumeratorHalfEdge {
                        edge: *edge_id,
                        flow: Flow::Source,
                    },
                    edge.data.pdg,
                    edge.directed.then_some(true),
                ));
            }
            if edge.vertices.1 == vertex {
                available.push((
                    NumeratorHalfEdge {
                        edge: *edge_id,
                        flow: Flow::Sink,
                    },
                    edge.data.pdg,
                    edge.directed.then_some(false),
                ));
            }
        }

        let mut ordered = Vec::with_capacity(rule.particles.len());
        for (leg, particle_name) in rule.particles.iter().enumerate() {
            let particle = self.model.particle(particle_name)?;
            let expected = ResolvedProcess::vertex_half_edge(&self.model, particle)?;
            let Some(position) = available.iter().position(|(_, pdg, direction)| {
                *pdg == expected.data.pdg && *direction == expected.direction
            }) else {
                return Err(GenerationError::MissingInteractionLeg {
                    vertex,
                    interaction: rule.name.clone(),
                    leg: leg + 1,
                    particle: particle_name.clone(),
                });
            };
            ordered.push(available.remove(position).0);
        }
        Ok(ordered)
    }

    fn propagator_numerator<'model>(
        &'model self,
        particle: &Particle,
    ) -> Result<&'model str, GenerationError> {
        if let Some(propagator) = &particle.propagator {
            return Ok(&self.model.propagator(propagator)?.numerator);
        }
        Ok(self
            .model
            .propagators()
            .iter()
            .find(|propagator| propagator.particle == particle.name)
            .map_or("1", |propagator| propagator.numerator.as_str()))
    }

    fn interaction_numerator(&self, rule: &VertexRule) -> Result<String, GenerationError> {
        let mut terms = Vec::new();
        for (color_index, color) in rule.color_structures.iter().enumerate() {
            for (lorentz_index, lorentz_name) in rule.lorentz_structures.iter().enumerate() {
                if let Some(coupling) = rule
                    .couplings
                    .get(color_index)
                    .and_then(|row| row.get(lorentz_index))
                    .and_then(Option::as_ref)
                {
                    let coupling = &self.model.coupling(coupling)?.expression;
                    let lorentz = &self.model.lorentz_structure(lorentz_name)?.structure;
                    terms.push(format!("({coupling})*({color})*({lorentz})"));
                }
            }
        }
        if terms.is_empty() {
            Ok("1".to_owned())
        } else {
            Ok(terms.join("+"))
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct EdgeColor {
    pdg: i64,
    direction: Option<bool>,
}

impl fmt::Display for EdgeColor {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.pdg.fmt(formatter)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct ExternalNode {
    index: usize,
    state: ExternalState,
    particle: ParticleReference,
    symmetry_class: usize,
}

#[derive(Debug, Clone, Default)]
struct NodeColor {
    external: Option<ExternalNode>,
}

impl NodeColor {
    fn key(&self) -> Option<(ExternalState, &ParticleReference, usize)> {
        self.external
            .as_ref()
            .map(|external| (external.state, &external.particle, external.symmetry_class))
    }
}

impl PartialEq for NodeColor {
    fn eq(&self, other: &Self) -> bool {
        self.key() == other.key()
    }
}

impl Eq for NodeColor {}

impl PartialOrd for NodeColor {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for NodeColor {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.key().cmp(&other.key())
    }
}

impl Hash for NodeColor {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.key().hash(state);
    }
}

impl fmt::Display for NodeColor {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.external {
            Some(external) => write!(formatter, "ext{}", external.index),
            None => formatter.write_str("v"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum ColoredNode {
    External(ExternalNode),
    Interaction(String),
}

type InteractionAssignment = (Graph<ColoredNode, EdgeColor>, usize);

impl fmt::Display for ColoredNode {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::External(external) => write!(formatter, "ext{}", external.index),
            Self::Interaction(rule) => formatter.write_str(rule),
        }
    }
}

struct ColoredTopology {
    graph: Graph<ColoredNode, EdgeColor>,
    symmetry: u64,
    multiplicity: usize,
    cut_partitions: Vec<CutPartition>,
}

struct CutPartition {
    left: CutSide,
    right: CutSide,
}

struct CutSide {
    coupling_orders: BTreeMap<String, usize>,
    loop_count: usize,
}

struct UnresolvedCutContent {
    particles: BTreeSet<i64>,
    maximum_multiplicity: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct FermionSigns {
    include_external_ordering: bool,
    external_ordering_negative: bool,
    include_antifermion_spin_sum: bool,
    antifermion_spin_sum_negative: bool,
}

struct NormalizedFermionFlow {
    graph: Graph<ColoredNode, EdgeColor>,
    signs: FermionSigns,
}

fn normalize_fermion_component(
    graph: &Graph<ColoredNode, EdgeColor>,
    starting_edge: usize,
    adjacency: &[Vec<(usize, usize)>],
    visited_edges: &mut [bool],
    normalized_vertices: &mut [Option<(usize, usize)>],
) -> BTreeSet<usize> {
    let (left, right) = graph.edges()[starting_edge].vertices;
    visited_edges[starting_edge] = true;
    normalized_vertices[starting_edge] = Some((left, right));
    let mut external_legs = BTreeSet::new();
    let mut record_external = |node_id: usize| {
        if let ColoredNode::External(external) = &graph.nodes()[node_id].data {
            external_legs.insert(external.index);
        }
    };
    record_external(left);
    record_external(right);

    // Read both ways from the first edge, flipping each subsequent edge when
    // needed so that the complete chain has one consistent orientation.
    for (mut current, read_forward) in [(right, true), (left, false)] {
        while let Some(&(edge_id, next)) = adjacency[current]
            .iter()
            .find(|(edge_id, _)| !visited_edges[*edge_id])
        {
            visited_edges[edge_id] = true;
            normalized_vertices[edge_id] = Some(if read_forward {
                (current, next)
            } else {
                (next, current)
            });
            record_external(next);
            current = next;
        }
    }
    external_legs
}

fn all_incoming_particle<'model>(
    external: &ExternalNode,
    model: &'model Model,
) -> Result<&'model Particle, GenerationError> {
    let particle = model.particle_by_pdg(external.particle.pdg)?;
    match external.state {
        ExternalState::Incoming => Ok(particle),
        ExternalState::Outgoing => Ok(model.antiparticle(particle)?),
    }
}

fn permutation_is_odd(values: &[usize]) -> bool {
    values
        .iter()
        .enumerate()
        .map(|(index, value)| {
            values[index + 1..]
                .iter()
                .filter(|next| value > *next)
                .count()
        })
        .sum::<usize>()
        % 2
        == 1
}

struct InteractionSignatures(BTreeMap<Vec<HalfEdge<EdgeColor>>, Vec<String>>);

impl InteractionSignatures {
    fn new(model: &Model, options: &GenerationOptions) -> Result<Self, GenerationError> {
        let allowed = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::VertexAllow(names) => Some(names.iter().collect::<BTreeSet<_>>()),
                _ => None,
            });
        let vetoed = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::VertexVeto(names) => Some(names.iter().collect::<BTreeSet<_>>()),
                _ => None,
            });
        let particle_vetoes = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::ParticleVeto(pdgs) => Some(pdgs.as_slice()),
                _ => None,
            });
        let mut signatures: BTreeMap<Vec<_>, Vec<_>> = BTreeMap::new();
        for rule in model.vertex_rules() {
            if allowed
                .as_ref()
                .is_some_and(|allowed| !allowed.contains(&rule.name))
                || vetoed
                    .as_ref()
                    .is_some_and(|vetoed| vetoed.contains(&rule.name))
            {
                continue;
            }
            let mut signature = Vec::new();
            let mut rejected = false;
            for particle_name in &rule.particles {
                let particle = model.particle(particle_name)?;
                if particle_vetoes.is_some_and(|vetoed| {
                    vetoed.contains(&particle.pdg_code)
                        || model
                            .antiparticle(particle)
                            .is_ok_and(|anti| vetoed.contains(&anti.pdg_code))
                }) {
                    rejected = true;
                    break;
                }
                signature.push(ResolvedProcess::vertex_half_edge(model, particle)?);
            }
            if !rejected {
                signature.sort();
                signatures
                    .entry(signature)
                    .or_default()
                    .push(rule.name.clone());
            }
        }
        // Stable signature and rule ordering makes graph generation
        // reproducible.
        for rules in signatures.values_mut() {
            rules.sort();
        }
        Ok(Self(signatures))
    }

    fn keys(&self) -> impl Iterator<Item = &Vec<HalfEdge<EdgeColor>>> {
        self.0.keys()
    }

    fn rules(&self, signature: &[EdgeColor]) -> Option<&[String]> {
        let mut half_edges: Vec<_> = signature
            .iter()
            .map(|edge| HalfEdge {
                direction: edge.direction,
                data: EdgeColor {
                    pdg: edge.pdg,
                    direction: None,
                },
            })
            .collect();
        // `EdgeColor` and `HalfEdge<EdgeColor>` have different derived orderings:
        // the former orders by PDG first, while the latter orders by direction
        // first. Restore the key ordering used when interaction signatures were
        // built before performing the lookup.
        half_edges.sort();
        self.0.get(&half_edges).map(Vec::as_slice)
    }
}

struct ResolvedProcess {
    generation_type: GenerationType,
    incoming: Vec<Particle>,
    outgoing: Vec<Vec<Particle>>,
    symmetrize_initial: bool,
    symmetrize_final: bool,
    symmetrize_left_right: bool,
    symmetrize_external_fermions: bool,
}

impl ResolvedProcess {
    fn new(model: &Model, process: &Process) -> Result<Self, GenerationError> {
        let resolve = |selector: &ParticleSelector| -> Result<Particle, ModelError> {
            match selector {
                ParticleSelector::Name(name) => model.particle(name).cloned(),
                ParticleSelector::Pdg(pdg) => model.particle_by_pdg(*pdg).cloned(),
            }
        };
        Ok(Self {
            generation_type: process.generation_type(),
            incoming: process
                .incoming()
                .iter()
                .map(resolve)
                .collect::<Result<_, _>>()?,
            outgoing: process
                .outgoing_alternatives()
                .iter()
                .map(|alternative| alternative.iter().map(resolve).collect())
                .collect::<Result<_, _>>()?,
            symmetrize_initial: process.symmetrizes_initial(),
            symmetrize_final: process.symmetrizes_final(),
            symmetrize_left_right: process.symmetrizes_left_right(),
            symmetrize_external_fermions: process.symmetrizes_external_fermions(),
        })
    }

    fn external_edges(
        &self,
        model: &Model,
    ) -> Result<Vec<(NodeColor, HalfEdge<EdgeColor>)>, GenerationError> {
        let mut result = Vec::new();
        for (index, particle) in self.incoming.iter().enumerate() {
            result.push(self.external_edge(model, particle, index, ExternalState::Incoming)?);
        }
        match self.generation_type {
            GenerationType::Amplitude => {
                for (offset, particle) in self.outgoing[0].iter().enumerate() {
                    result.push(self.external_edge(
                        model,
                        particle,
                        self.incoming.len() + offset,
                        ExternalState::Outgoing,
                    )?);
                }
            }
            GenerationType::CrossSection => {
                for (offset, particle) in self.incoming.iter().enumerate() {
                    result.push(self.external_edge(
                        model,
                        particle,
                        self.incoming.len() + offset,
                        ExternalState::Outgoing,
                    )?);
                }
            }
        }
        Ok(result)
    }

    fn external_edge(
        &self,
        model: &Model,
        particle: &Particle,
        index: usize,
        state: ExternalState,
    ) -> Result<(NodeColor, HalfEdge<EdgeColor>), GenerationError> {
        // Symbolica represents every non-self-conjugate species by its particle
        // entry; particle versus antiparticle is encoded by edge direction.
        let base = if particle.is_antiparticle() {
            model.antiparticle(particle)?
        } else {
            particle
        };
        let direction = if particle.is_self_antiparticle() {
            None
        } else {
            Some(match (state, particle.is_antiparticle()) {
                (ExternalState::Incoming, false) | (ExternalState::Outgoing, true) => true,
                (ExternalState::Incoming, true) | (ExternalState::Outgoing, false) => false,
            })
        };
        let symmetrized = (match state {
            ExternalState::Incoming => self.symmetrize_initial,
            ExternalState::Outgoing => self.symmetrize_final,
        } || self.symmetrize_left_right)
            && (self.generation_type == GenerationType::CrossSection
                || !particle.is_fermion()
                || self.symmetrize_external_fermions);
        Ok((
            NodeColor {
                external: Some(ExternalNode {
                    index,
                    state,
                    particle: ParticleReference::new(&particle.name, particle.pdg_code),
                    symmetry_class: if symmetrized { 0 } else { index + 1 },
                }),
            },
            HalfEdge {
                direction,
                data: EdgeColor {
                    pdg: base.pdg_code,
                    direction: None,
                },
            },
        ))
    }

    // Normalize every fermion chain before exposing the diagram. A chain keeps
    // the orientation of its first external edge and all following edges are
    // oriented consistently with it. Purely virtual loops start from their
    // first edge in canonical edge order.
    fn normalize_fermion_flows(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
        model: &Model,
    ) -> Result<NormalizedFermionFlow, GenerationError> {
        let mut adjacency = vec![Vec::new(); graph.nodes().len()];
        let mut degrees = vec![0; graph.nodes().len()];
        for (edge_id, edge) in graph.edges().iter().enumerate() {
            if !model.particle_by_pdg(edge.data.pdg)?.is_fermion() {
                continue;
            }
            let (left, right) = edge.vertices;
            adjacency[left].push((edge_id, right));
            if left == right {
                degrees[left] += 2;
            } else {
                adjacency[right].push((edge_id, left));
                degrees[left] += 1;
                degrees[right] += 1;
            }
        }
        for (vertex, degree) in degrees.iter().copied().enumerate() {
            if degree > 2 {
                return Err(GenerationError::UnsupportedFermionBranching { vertex, degree });
            }
        }

        let external_nodes: BTreeMap<_, _> = graph
            .nodes()
            .iter()
            .enumerate()
            .filter_map(|(node_id, node)| match &node.data {
                ColoredNode::External(external) => {
                    Some((external.index, (node_id, external.clone())))
                }
                ColoredNode::Interaction(_) => None,
            })
            .collect();
        let mut visited_edges = vec![false; graph.edges().len()];
        let mut normalized_vertices = vec![None; graph.edges().len()];
        // Pairing of the external fermion flows. The keys are sorted pairs of
        // all-incoming PDGs, and values are connected external leg IDs in the
        // order of that key.
        let mut pairings: BTreeMap<(i64, i64), Vec<(usize, usize)>> = BTreeMap::new();
        let mut line_pairings = Vec::new();

        // First fix flows connected to external legs, then all internal flows.
        for (node_id, _) in external_nodes.values() {
            let Some((edge_id, _)) = adjacency[*node_id]
                .iter()
                .find(|(edge_id, _)| !visited_edges[*edge_id])
            else {
                continue;
            };
            let legs = normalize_fermion_component(
                graph,
                *edge_id,
                &adjacency,
                &mut visited_edges,
                &mut normalized_vertices,
            );
            let legs = legs.into_iter().collect::<Vec<_>>();
            let [first_leg, second_leg] = legs.as_slice() else {
                return Err(GenerationError::InvalidExternalFermionLegCount { legs });
            };
            let first = external_nodes
                .get(first_leg)
                .ok_or(GenerationError::MissingExternalParticle { index: *first_leg })?;
            let second = external_nodes
                .get(second_leg)
                .ok_or(GenerationError::MissingExternalParticle { index: *second_leg })?;
            let first_particle = all_incoming_particle(&first.1, model)?;
            let second_particle = all_incoming_particle(&second.1, model)?;
            let (anti_leg, anti, particle_leg, particle) = match (
                first_particle.is_fermion(),
                first_particle.is_antiparticle(),
                second_particle.is_fermion(),
                second_particle.is_antiparticle(),
            ) {
                (true, true, true, false) => {
                    (*first_leg, first_particle, *second_leg, second_particle)
                }
                (true, false, true, true) => {
                    (*second_leg, second_particle, *first_leg, first_particle)
                }
                _ => {
                    return Err(GenerationError::InvalidExternalFermionPair {
                        legs: [*first_leg, *second_leg],
                        particles: [first_particle.pdg_code, second_particle.pdg_code],
                    });
                }
            };
            pairings
                .entry((anti.pdg_code, particle.pdg_code))
                .or_default()
                .push((anti_leg, particle_leg));
            line_pairings.push((*first_leg, *second_leg));
        }
        for edge_id in 0..graph.edges().len() {
            if model
                .particle_by_pdg(graph.edges()[edge_id].data.pdg)?
                .is_fermion()
                && !visited_edges[edge_id]
            {
                normalize_fermion_component(
                    graph,
                    edge_id,
                    &adjacency,
                    &mut visited_edges,
                    &mut normalized_vertices,
                );
            }
        }

        let mut normalized_graph = Graph::new();
        for node in graph.nodes() {
            normalized_graph.add_node(node.data.clone());
        }
        for (edge_id, edge) in graph.edges().iter().enumerate() {
            let (left, right) = normalized_vertices[edge_id].unwrap_or(edge.vertices);
            normalized_graph
                .add_edge(left, right, edge.directed, edge.data)
                .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
        }

        let concatenated_lines = pairings
            .values()
            .flat_map(|lines| lines.iter().flat_map(|(anti, particle)| [*anti, *particle]))
            .collect::<Vec<_>>();
        let external_ordering_negative = match self.generation_type {
            GenerationType::Amplitude
                if self.symmetrize_external_fermions
                    && (self.symmetrize_initial || self.symmetrize_final) =>
            {
                false
            }
            GenerationType::Amplitude => permutation_is_odd(&concatenated_lines),
            GenerationType::CrossSection => {
                self.external_fermion_loop_count(&line_pairings)? % 2 == 1
            }
        };
        let antifermion_spin_sum_negative = self.generation_type == GenerationType::CrossSection
            && self
                .incoming
                .iter()
                .filter(|particle| particle.is_fermion() && particle.is_antiparticle())
                .count()
                % 2
                == 1;

        Ok(NormalizedFermionFlow {
            graph: normalized_graph,
            signs: FermionSigns {
                include_external_ordering: self.generation_type == GenerationType::CrossSection
                    || !(self.symmetrize_external_fermions
                        && (self.symmetrize_initial || self.symmetrize_final)),
                external_ordering_negative,
                include_antifermion_spin_sum: self.generation_type == GenerationType::CrossSection,
                antifermion_spin_sum_negative,
            },
        })
    }

    fn external_fermion_loop_count(
        &self,
        line_pairings: &[(usize, usize)],
    ) -> Result<usize, GenerationError> {
        let mut adjacency: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        let mut connect = |left: usize, right: usize| {
            adjacency.entry(left).or_default().push(right);
            adjacency.entry(right).or_default().push(left);
        };
        for &(left, right) in line_pairings {
            connect(left, right);
        }
        for (index, particle) in self.incoming.iter().enumerate() {
            if particle.is_fermion() {
                connect(index, self.incoming.len() + index);
            }
        }
        for (&leg, neighbors) in &adjacency {
            if neighbors.len() != 2 {
                return Err(GenerationError::InvalidExternalFermionLoopDegree {
                    leg,
                    degree: neighbors.len(),
                });
            }
        }

        let mut visited = BTreeSet::new();
        let mut loops = 0;
        for &start in adjacency.keys() {
            if !visited.insert(start) {
                continue;
            }
            loops += 1;
            let mut stack = vec![start];
            while let Some(leg) = stack.pop() {
                if let Some(neighbors) = adjacency.get(&leg) {
                    for &neighbor in neighbors {
                        if visited.insert(neighbor) {
                            stack.push(neighbor);
                        }
                    }
                }
            }
        }
        Ok(loops)
    }

    fn vertex_half_edge(
        model: &Model,
        particle: &Particle,
    ) -> Result<HalfEdge<EdgeColor>, GenerationError> {
        let base = if particle.is_antiparticle() {
            model.antiparticle(particle)?
        } else {
            particle
        };
        Ok(HalfEdge {
            direction: if particle.is_self_antiparticle() {
                None
            } else {
                Some(particle.is_antiparticle())
            },
            data: EdgeColor {
                pdg: base.pdg_code,
                direction: None,
            },
        })
    }

    // Enumerate exact s-channel cuts for cross-section generation. The caller
    // applies amplitude-level constraints to the returned partitions.
    fn cut_partitions(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
        model: &Model,
        options: &GenerationOptions,
    ) -> Result<Vec<CutPartition>, GenerationError> {
        let unresolved = Self::unresolved_cut_content(model, options)?;
        let targets: Vec<Vec<i64>> = self
            .outgoing
            .iter()
            .map(|particles| {
                let mut pdgs = particles
                    .iter()
                    .map(|particle| particle.pdg_code)
                    .collect::<Vec<_>>();
                pdgs.sort();
                pdgs
            })
            .collect();
        let he_graph: HedgeGraph<EdgeColor, ColoredNode, ()> = graph.clone().into();
        let mut source = Vec::new();
        let mut target = Vec::new();
        for (node_id, _, node) in he_graph.iter_nodes() {
            if let ColoredNode::External(external) = node {
                match external.state {
                    ExternalState::Incoming => source.push(node_id),
                    ExternalState::Outgoing => target.push(node_id),
                }
            }
        }
        if source.is_empty() || target.is_empty() {
            // Vacuum or degenerate graphs have no physical s-channel cut.
            return Ok(Vec::new());
        }

        let default_blob_range = 1..=1;
        let default_spectator_range = 0..=0;
        let blob_range = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::BlobRange(range) => Some(range),
                _ => None,
            })
            .unwrap_or(&default_blob_range);
        let spectator_range = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::SpectatorRange(range) => Some(range),
                _ => None,
            })
            .unwrap_or(&default_spectator_range);

        let validate_connectivity = |subgraph: &SuBitGraph| {
            let mut blobs = 0;
            let mut spectators = 0;
            for component in he_graph.connected_components(subgraph) {
                // A component with more than one included half-edge is a blob;
                // a singleton component is a spectator.
                if component.n_included() > 1 {
                    blobs += 1;
                } else {
                    spectators += 1;
                }
            }
            blob_range.contains(&blobs) && spectator_range.contains(&spectators)
        };

        let summarize_side = |subgraph: &SuBitGraph| -> Result<CutSide, GenerationError> {
            let mut coupling_orders = BTreeMap::new();
            for (_, _, node) in he_graph.iter_nodes_of(subgraph) {
                let ColoredNode::Interaction(rule_name) = node else {
                    continue;
                };
                for (name, order) in model.vertex_rule(rule_name)?.coupling_orders(model)? {
                    *coupling_orders.entry(name).or_insert(0) += order;
                }
            }
            let internal = InternalSubGraph::cleaned_filter_pessimist(subgraph.clone(), &he_graph);
            Ok(CutSide {
                coupling_orders,
                loop_count: he_graph.cyclotomatic_number(&internal),
            })
        };

        let mut partitions = Vec::new();
        for (left, cut, right) in he_graph.all_cuts_from_ids(&source, &target) {
            if !validate_connectivity(&left) || !validate_connectivity(&right) {
                continue;
            }
            // Convert sink-oriented cut half-edges to antiparticles so the
            // content describes the physical state crossing left to right.
            let mut cut_content = cut
                .iter_left_hedges()
                .map(|hedge| {
                    let particle = model.particle_by_pdg(he_graph[[&hedge]].pdg)?;
                    if matches!(he_graph.flow(hedge), Flow::Sink) {
                        model.antiparticle(particle).map(|anti| anti.pdg_code)
                    } else {
                        Ok(particle.pdg_code)
                    }
                })
                .collect::<Result<Vec<_>, ModelError>>()?;
            cut_content.sort();
            if !targets
                .iter()
                .any(|target| Self::matches_cut_content(&cut_content, target, unresolved.as_ref()))
            {
                continue;
            }
            partitions.push(CutPartition {
                left: summarize_side(&left)?,
                right: summarize_side(&right)?,
            });
        }
        Ok(partitions)
    }

    fn unresolved_cut_content(
        model: &Model,
        options: &GenerationOptions,
    ) -> Result<Option<UnresolvedCutContent>, GenerationError> {
        let Some(required) = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::PerturbativeOrders(required) => Some(required),
                _ => None,
            })
        else {
            return Ok(None);
        };
        let mut particles = BTreeSet::new();
        for rule in model.vertex_rules() {
            let orders = rule.coupling_orders(model)?;
            if !required.keys().any(|name| orders.contains_key(name)) {
                continue;
            }
            for particle_name in &rule.particles {
                let particle = model.particle(particle_name)?;
                if particle.is_massless() {
                    particles.insert(particle.pdg_code);
                }
            }
        }
        let maximum_multiplicity = required.values().try_fold(0_usize, |sum, order| {
            sum.checked_add(*order)
                .ok_or(GenerationError::ArithmeticOverflow(
                    "the unresolved cut multiplicity",
                ))
        })?;
        Ok(Some(UnresolvedCutContent {
            particles,
            maximum_multiplicity,
        }))
    }

    fn matches_cut_content(
        cut_content: &[i64],
        target: &[i64],
        unresolved: Option<&UnresolvedCutContent>,
    ) -> bool {
        let mut remainder = cut_content.to_vec();
        for particle in target {
            let Some(position) = remainder.iter().position(|candidate| candidate == particle)
            else {
                return false;
            };
            remainder.remove(position);
        }
        remainder.is_empty()
            || unresolved.is_some_and(|unresolved| {
                remainder.len() <= unresolved.maximum_multiplicity
                    && remainder
                        .iter()
                        .all(|particle| unresolved.particles.contains(particle))
            })
    }

    fn passes_sewn_filter(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
        options: &GenerationOptions,
    ) -> bool {
        let Some(SewnFilterOptions {
            filter_tadpoles: true,
        }) = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::Sewn(options) => Some(*options),
                _ => None,
            })
        else {
            return true;
        };

        let mut graph: HedgeGraph<EdgeColor, ColoredNode, ()> = graph.clone().into();
        let external_nodes: Vec<_> = graph
            .iter_nodes()
            .filter_map(|(node, _, color)| match color {
                ColoredNode::External(external) => Some((external.index, node, color.clone())),
                ColoredNode::Interaction(_) => None,
            })
            .collect();
        let externals: Vec<_> = external_nodes.iter().map(|(_, node, _)| *node).collect();
        if externals.is_empty() {
            return true;
        }
        let connected_components_before = graph.tadpoles(&externals).len() + 1;
        let by_index: BTreeMap<usize, (NodeIndex, ColoredNode)> = external_nodes
            .into_iter()
            .map(|(index, node, color)| (index, (node, color)))
            .collect();
        for incoming in 0..self.incoming.len() {
            let outgoing = self.incoming.len() + incoming;
            if let (Some((left, color)), Some((right, _))) =
                (by_index.get(&incoming), by_index.get(&outgoing))
            {
                graph.identify_nodes(&[*left, *right], color.clone());
            }
        }
        // NodeStorageVec retains the crowns of nodes that were identified so
        // callers can inspect the identification history. Cycle detection must
        // run on the sewn quotient graph, without those historical crowns.
        graph.forget_identification_history();
        let non_bridges = graph.non_bridges();
        graph.count_connected_components(&non_bridges) == connected_components_before
    }
}

#[cfg(test)]
mod tests {
    use feynkit_model::{
        ComplexValue, Coupling, LorentzStructure, ModelDefinition, Order, Parameter,
        ParameterNature, ParameterType, Propagator,
    };

    use crate::{CancellationToken, GraphGroupingOptions, NumeratorGrouping};

    use super::*;

    fn scalar_model() -> Model {
        Model::new(ModelDefinition {
            name: "phi3".to_owned(),
            restriction: None,
            orders: vec![Order {
                name: "QED".to_owned(),
                expansion_order: 99,
                hierarchy: 1,
            }],
            parameters: vec![
                Parameter {
                    name: "ZERO".to_owned(),
                    lhablock: None,
                    lhacode: None,
                    nature: ParameterNature::Internal,
                    parameter_type: ParameterType::Real,
                    value: Some(ComplexValue::ZERO),
                    expression: None,
                },
                Parameter {
                    name: "M".to_owned(),
                    lhablock: Some("MASS".to_owned()),
                    lhacode: Some(vec![25]),
                    nature: ParameterNature::External,
                    parameter_type: ParameterType::Real,
                    value: Some(ComplexValue::new(1.0, 0.0)),
                    expression: None,
                },
            ],
            particles: vec![Particle {
                pdg_code: 25,
                name: "phi".to_owned(),
                antiname: "phi".to_owned(),
                spin: 1,
                color: 1,
                mass: "M".to_owned(),
                width: "ZERO".to_owned(),
                texname: "phi".to_owned(),
                antitexname: "phi".to_owned(),
                charge: 0.0,
                ghost_number: 0,
                lepton_number: 0,
                y_charge: 0,
                propagating: true,
                goldstone: false,
                propagator: None,
            }],
            propagators: vec![Propagator {
                name: "phi_prop".to_owned(),
                particle: "phi".to_owned(),
                numerator: "1".to_owned(),
                denominator: "P^2-M^2".to_owned(),
            }],
            lorentz_structures: vec![LorentzStructure {
                name: "L1".to_owned(),
                spins: vec![1, 1, 1],
                structure: "1".to_owned(),
            }],
            couplings: vec![Coupling {
                name: "GC1".to_owned(),
                expression: "g".to_owned(),
                orders: BTreeMap::from([("QED".to_owned(), 1)]),
                value: None,
            }],
            vertex_rules: vec![VertexRule {
                name: "V1".to_owned(),
                particles: vec!["phi".to_owned(); 3],
                color_structures: vec!["1".to_owned()],
                lorentz_structures: vec!["L1".to_owned()],
                couplings: vec![vec![Some("GC1".to_owned())]],
            }],
            functions: Vec::new(),
            form_factors: Vec::new(),
        })
        .unwrap()
    }

    fn fermion_model() -> Model {
        Model::from_json(
            r#"{
                "name":"toy_qed",
                "restriction":null,
                "orders":[
                    {"name":"QED","expansion_order":99,"hierarchy":1},
                    {"name":"QCD","expansion_order":99,"hierarchy":1}
                ],
                "parameters":[{
                    "name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal",
                    "parameter_type":"real","value":[0.0,0.0],"expression":null
                }],
                "particles":[
                    {"pdg_code":1,"name":"f","antiname":"f~","spin":2,"color":1,
                     "mass":"ZERO","width":"ZERO","texname":"f","antitexname":"f~",
                     "charge":-1.0,"ghost_number":0,"lepton_number":1,"y_charge":0},
                    {"pdg_code":-1,"name":"f~","antiname":"f","spin":2,"color":1,
                     "mass":"ZERO","width":"ZERO","texname":"f~","antitexname":"f",
                     "charge":1.0,"ghost_number":0,"lepton_number":-1,"y_charge":0},
                    {"pdg_code":22,"name":"a","antiname":"a","spin":3,"color":1,
                     "mass":"ZERO","width":"ZERO","texname":"a","antitexname":"a",
                     "charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0}
                ],
                "propagators":[
                    {"name":"f_prop","particle":"f","numerator":"Slash(P)","denominator":"P^2"},
                    {"name":"a_prop","particle":"a","numerator":"1","denominator":"P^2"}
                ],
                "lorentz_structures":[{
                    "name":"FFV","spins":[2,2,3],"structure":"Gamma(3,2,1)"
                }],
                "couplings":[
                    {"name":"GC_LOW","expression":"g_low","orders":[["QED",1]],"value":null},
                    {"name":"GC_HIGH","expression":"g_high","orders":[["QED",3],["QCD",2]],"value":null}
                ],
                "vertex_rules":[{
                    "name":"V","particles":["f~","f","a"],
                    "color_structures":["1","T(1,2,3)"],
                    "lorentz_structures":["FFV"],
                    "couplings":[["GC_LOW"],["GC_HIGH"]]
                }]
            }"#,
        )
        .unwrap()
    }

    fn exact_cut_model() -> Model {
        Model::from_json(
            r#"{
                "name":"exact_cut_regression",
                "restriction":null,
                "orders":[
                    {"name":"QED","expansion_order":99,"hierarchy":2},
                    {"name":"QCD","expansion_order":99,"hierarchy":1}
                ],
                "parameters":[
                    {"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal",
                     "parameter_type":"real","value":[0.0,0.0],"expression":null},
                    {"name":"M","lhablock":"MASS","lhacode":[1],"nature":"external",
                     "parameter_type":"real","value":[1.0,0.0],"expression":null}
                ],
                "particles":[
                    {"pdg_code":22,"name":"a","antiname":"a","spin":3,"color":1,
                     "mass":"ZERO","width":"ZERO","texname":"a","antitexname":"a",
                     "charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":21,"name":"g","antiname":"g","spin":3,"color":8,
                     "mass":"ZERO","width":"ZERO","texname":"g","antitexname":"g",
                     "charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":5,"name":"b","antiname":"b~","spin":2,"color":3,
                     "mass":"M","width":"ZERO","texname":"b","antitexname":"b~",
                     "charge":-0.3333333333333333,"ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":-5,"name":"b~","antiname":"b","spin":2,"color":-3,
                     "mass":"M","width":"ZERO","texname":"b~","antitexname":"b",
                     "charge":0.3333333333333333,"ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":25,"name":"H","antiname":"H","spin":1,"color":1,
                     "mass":"M","width":"ZERO","texname":"H","antitexname":"H",
                     "charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":11,"name":"e-","antiname":"e+","spin":2,"color":1,
                     "mass":"M","width":"ZERO","texname":"e-","antitexname":"e+",
                     "charge":-1.0,"ghost_number":0,"lepton_number":1,"y_charge":0},
                    {"pdg_code":-11,"name":"e+","antiname":"e-","spin":2,"color":1,
                     "mass":"M","width":"ZERO","texname":"e+","antitexname":"e-",
                     "charge":1.0,"ghost_number":0,"lepton_number":-1,"y_charge":0}
                ],
                "propagators":[],
                "lorentz_structures":[
                    {"name":"FFV","spins":[2,2,3],"structure":"Gamma(3,2,1)"},
                    {"name":"FFS","spins":[2,2,1],"structure":"1"}
                ],
                "couplings":[
                    {"name":"GC_QED","expression":"ee","orders":[["QED",1]],"value":null},
                    {"name":"GC_QCD","expression":"gs","orders":[["QCD",1]],"value":null}
                ],
                "vertex_rules":[
                    {"name":"V_73","particles":["b~","b","a"],"color_structures":["1"],
                     "lorentz_structures":["FFV"],"couplings":[["GC_QED"]]},
                    {"name":"V_76","particles":["b~","b","g"],"color_structures":["1"],
                     "lorentz_structures":["FFV"],"couplings":[["GC_QCD"]]},
                    {"name":"V_78","particles":["b~","b","H"],"color_structures":["1"],
                     "lorentz_structures":["FFS"],"couplings":[["GC_QED"]]},
                    {"name":"V_98","particles":["e+","e-","a"],"color_structures":["1"],
                     "lorentz_structures":["FFV"],"couplings":[["GC_QED"]]}
                ]
            }"#,
        )
        .unwrap()
    }

    fn photon_process() -> Process {
        Process::amplitude([22_i64], [22_i64])
    }

    #[test]
    fn generates_tree_and_one_loop_diagrams() {
        let generator = Generator::new(scalar_model());
        let process = Process::amplitude(["phi"], ["phi", "phi"])
            .with_loop_count(0, 1)
            .unwrap();
        let generated = generator
            .generate(
                &process,
                &GenerationOptions::default().allow_self_loops(true),
            )
            .unwrap();
        assert!(!generated.diagrams.is_empty());
        assert!(generated.report.completed);
        assert!(
            generated
                .diagrams
                .iter()
                .all(|diagram| diagram.numerator().is_some())
        );
    }

    #[test]
    fn cancelled_generation_returns_an_incomplete_result() {
        let cancellation = CancellationToken::new();
        cancellation.cancel();
        let generated = Generator::new(scalar_model())
            .generate(
                &Process::amplitude(["phi"], ["phi", "phi"]),
                &GenerationOptions::default().cancellation_token(cancellation),
            )
            .unwrap();

        assert!(!generated.report.completed);
        assert!(generated.diagrams.is_empty());
    }

    #[test]
    fn cancellation_callback_returns_an_incomplete_result() {
        let generated = Generator::new(scalar_model())
            .generate(
                &Process::amplitude(["phi"], ["phi", "phi"]),
                &GenerationOptions::default().cancellation_check(|| true),
            )
            .unwrap();

        assert!(!generated.report.completed);
        assert!(generated.diagrams.is_empty());
    }

    #[test]
    fn directed_and_undirected_half_edges_use_the_interaction_key_order() {
        let model = fermion_model();
        let options = GenerationOptions::default()
            .with_graph_filter(GenerationFilter::VertexAllow(vec!["V".to_owned()]));
        let signatures = InteractionSignatures::new(&model, &options).unwrap();
        let mut signature = vec![
            EdgeColor {
                pdg: 1,
                direction: Some(false),
            },
            EdgeColor {
                pdg: 1,
                direction: Some(true),
            },
            EdgeColor {
                pdg: 22,
                direction: None,
            },
        ];
        signature.sort();

        assert_eq!(signatures.rules(&signature).unwrap(), &["V".to_owned()]);
    }

    #[test]
    fn external_fermion_chains_are_oriented_consistently() {
        let model = fermion_model();
        let process =
            ResolvedProcess::new(&model, &Process::amplitude([1_i64, -1], Vec::<i64>::new()))
                .unwrap();
        let mut graph = Graph::new();
        let particle = graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Incoming,
            particle: ParticleReference::new("f", 1),
            symmetry_class: 1,
        }));
        let left = graph.add_node(ColoredNode::Interaction("V".to_owned()));
        let right = graph.add_node(ColoredNode::Interaction("V".to_owned()));
        let antiparticle = graph.add_node(ColoredNode::External(ExternalNode {
            index: 1,
            state: ExternalState::Incoming,
            particle: ParticleReference::new("f~", -1),
            symmetry_class: 2,
        }));
        let fermion = EdgeColor {
            pdg: 1,
            direction: None,
        };
        graph.add_edge(particle, left, true, fermion).unwrap();
        graph.add_edge(right, left, true, fermion).unwrap();
        graph.add_edge(right, antiparticle, true, fermion).unwrap();

        let normalized = process.normalize_fermion_flows(&graph, &model).unwrap();
        assert_eq!(normalized.graph.edges()[0].vertices, (particle, left));
        assert_eq!(normalized.graph.edges()[1].vertices, (left, right));
        assert_eq!(normalized.graph.edges()[2].vertices, (right, antiparticle));
        assert!(normalized.signs.external_ordering_negative);
        assert!(normalized.signs.include_external_ordering);
        assert!(!normalized.signs.include_antifermion_spin_sum);
        assert!(!normalized.signs.antifermion_spin_sum_negative);

        let symmetrized = ResolvedProcess::new(
            &model,
            &Process::amplitude([1_i64, -1], Vec::<i64>::new())
                .symmetrize_initial(true)
                .symmetrize_external_fermions(true),
        )
        .unwrap()
        .normalize_fermion_flows(&graph, &model)
        .unwrap();
        assert!(!symmetrized.signs.include_external_ordering);
        assert!(!symmetrized.signs.external_ordering_negative);
    }

    #[test]
    fn vertex_orders_take_each_order_maximum_and_graph_orders_sum_vertices() {
        let model = fermion_model();
        let rule_orders = model
            .vertex_rule("V")
            .unwrap()
            .coupling_orders(&model)
            .unwrap();
        assert_eq!(
            rule_orders,
            BTreeMap::from([("QCD".to_owned(), 2), ("QED".to_owned(), 3)])
        );

        let generator = Generator::new(model);
        let mut graph = Graph::new();
        graph.add_node(ColoredNode::Interaction("V".to_owned()));
        graph.add_node(ColoredNode::Interaction("V".to_owned()));
        assert_eq!(
            generator.coupling_orders(&graph).unwrap(),
            BTreeMap::from([("QCD".to_owned(), 4), ("QED".to_owned(), 6)])
        );
    }

    #[test]
    fn interaction_numerator_uses_model_expressions() {
        let model = fermion_model();
        let generator = Generator::new(model.clone());
        assert_eq!(
            generator
                .interaction_numerator(model.vertex_rule("V").unwrap())
                .unwrap(),
            "(g_low)*(1)*(Gamma(3,2,1))+(g_high)*(T(1,2,3))*(Gamma(3,2,1))"
        );
    }

    #[test]
    fn closed_fermion_loops_always_contribute_their_sign() {
        let generator = Generator::new(fermion_model());
        let mut graph = Graph::new();
        let interaction = graph.add_node(ColoredNode::Interaction("V".to_owned()));
        let external = graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Outgoing,
            particle: ParticleReference::new("a", 22),
            symmetry_class: 1,
        }));
        graph
            .add_edge(
                interaction,
                interaction,
                true,
                EdgeColor {
                    pdg: 1,
                    direction: None,
                },
            )
            .unwrap();
        graph
            .add_edge(
                interaction,
                external,
                false,
                EdgeColor {
                    pdg: 22,
                    direction: None,
                },
            )
            .unwrap();

        let loop_count = generator.closed_fermion_loop_count(&graph).unwrap();
        assert_eq!(loop_count, 1);
        let diagram = generator
            .to_diagram(
                "loop".to_owned(),
                ColoredTopology {
                    graph,
                    symmetry: 3,
                    multiplicity: 2,
                    cut_partitions: Vec::new(),
                },
                loop_count,
                FermionSigns::default(),
            )
            .unwrap();
        assert_eq!(
            diagram.overall_factor(),
            "InternalFermionLoopSign(-1)*CouplingsMultiplicity(2)/AutG(3)"
        );
    }

    #[test]
    fn positive_external_fermion_signs_keep_symbolic_provenance() {
        let generator = Generator::new(fermion_model());
        let diagram = generator
            .to_diagram(
                "positive".to_owned(),
                ColoredTopology {
                    graph: Graph::new(),
                    symmetry: 1,
                    multiplicity: 1,
                    cut_partitions: Vec::new(),
                },
                0,
                FermionSigns {
                    include_external_ordering: true,
                    external_ordering_negative: false,
                    include_antifermion_spin_sum: true,
                    antifermion_spin_sum_negative: false,
                },
            )
            .unwrap();
        assert_eq!(
            diagram.overall_factor(),
            "AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*1/AutG(1)"
        );
    }

    #[test]
    fn branching_fermion_chains_are_rejected() {
        let generator = Generator::new(fermion_model());
        let mut graph = Graph::new();
        let center = graph.add_node(ColoredNode::Interaction("V".to_owned()));
        for _ in 0..3 {
            let endpoint = graph.add_node(ColoredNode::Interaction("V".to_owned()));
            graph
                .add_edge(
                    center,
                    endpoint,
                    true,
                    EdgeColor {
                        pdg: 1,
                        direction: None,
                    },
                )
                .unwrap();
        }

        assert!(matches!(
            generator.closed_fermion_loop_count(&graph),
            Err(GenerationError::UnsupportedFermionBranching {
                vertex,
                degree: 3
            }) if vertex == center
        ));
    }

    #[test]
    fn filters_are_scoped_and_validated_before_generation() {
        let generator = Generator::new(fermion_model());
        let process = photon_process();

        let unsupported = GenerationOptions::default().with_graph_filter(
            GenerationFilter::SelfEnergy(crate::SelfEnergyFilterOptions {
                only_scaleless: true,
                ..crate::SelfEnergyFilterOptions::default()
            }),
        );
        assert!(matches!(
            generator.generate(&process, &unsupported),
            Err(GenerationError::UnsupportedFilterOption {
                scope: FilterScope::Graph,
                filter: GenerationFilterKind::SelfEnergy,
                option: "only_scaleless"
            })
        ));

        let perturbative_filter =
            GenerationFilter::PerturbativeOrders(BTreeMap::from([("QED".to_owned(), 1)]));
        generator
            .validate_filters(
                FilterScope::Graph,
                GenerationType::CrossSection,
                std::slice::from_ref(&perturbative_filter),
            )
            .unwrap();
        let perturbative = GenerationOptions::default().with_graph_filter(perturbative_filter);
        assert!(matches!(
            generator.generate(&process, &perturbative),
            Err(GenerationError::InvalidFilterScope {
                scope: FilterScope::Graph,
                filter: GenerationFilterKind::PerturbativeOrders,
                generation_type: GenerationType::Amplitude
            })
        ));

        let sewn_filter = GenerationFilter::Sewn(SewnFilterOptions::default());
        generator
            .validate_filters(
                FilterScope::Graph,
                GenerationType::CrossSection,
                std::slice::from_ref(&sewn_filter),
            )
            .unwrap();
        assert!(matches!(
            generator.generate(
                &process,
                &GenerationOptions::default().with_graph_filter(sewn_filter)
            ),
            Err(GenerationError::InvalidFilterScope {
                scope: FilterScope::Graph,
                filter: GenerationFilterKind::Sewn,
                generation_type: GenerationType::Amplitude
            })
        ));

        let duplicate = GenerationOptions::default()
            .with_graph_filter(GenerationFilter::LoopCountRange((0, 0)))
            .with_graph_filter(GenerationFilter::LoopCountRange((0, 1)));
        assert!(matches!(
            generator.generate(&process, &duplicate),
            Err(GenerationError::DuplicateFilter {
                scope: FilterScope::Graph,
                filter: GenerationFilterKind::LoopCountRange
            })
        ));

        let invalid = GenerationOptions::default()
            .with_graph_filter(GenerationFilter::LoopCountRange((2, 1)));
        assert!(matches!(
            generator.generate(&process, &invalid),
            Err(GenerationError::InvalidFilterRange {
                scope: FilterScope::Graph,
                filter: GenerationFilterKind::LoopCountRange,
                minimum: 2,
                maximum: 1
            })
        ));

        let invalid_scope = GenerationOptions::default()
            .with_cut_amplitude_filter(GenerationFilter::ParticleVeto(vec![1]));
        assert!(matches!(
            generator.generate(&process, &invalid_scope),
            Err(GenerationError::InvalidFilterScope {
                scope: FilterScope::CutAmplitude,
                filter: GenerationFilterKind::ParticleVeto,
                generation_type: GenerationType::Amplitude
            })
        ));
    }

    #[test]
    fn grouping_options_are_validated() {
        let generator = Generator::new(fermion_model());
        let process = photon_process();

        let grouping =
            GenerationOptions::default().numerator_grouping(NumeratorGrouping::OnlyDetectZeroes);
        generator.validate_options(&process, &grouping).unwrap();

        let invalid = GenerationOptions::default().numerator_grouping(
            NumeratorGrouping::Identical(GraphGroupingOptions {
                number_of_numerical_samples: 0,
                ..GraphGroupingOptions::default()
            }),
        );
        assert!(matches!(
            generator.validate_options(&process, &invalid),
            Err(GenerationError::InvalidNumericalSampleCount)
        ));
    }

    #[test]
    fn public_limits_report_arithmetic_overflow() {
        let model = fermion_model();
        let generator = Generator::new(model.clone());
        assert!(matches!(
            generator.generate(
                &photon_process(),
                &GenerationOptions::default().max_vertices(usize::MAX),
            ),
            Err(GenerationError::ArithmeticOverflow(
                "the Symbolica vertex limit"
            ))
        ));

        let options =
            GenerationOptions::default().with_graph_filter(GenerationFilter::PerturbativeOrders(
                BTreeMap::from([("QED".to_owned(), usize::MAX), ("QCD".to_owned(), 1)]),
            ));
        assert!(matches!(
            ResolvedProcess::unresolved_cut_content(&model, &options),
            Err(GenerationError::ArithmeticOverflow(
                "the unresolved cut multiplicity"
            ))
        ));
    }

    #[test]
    fn perturbative_orders_allow_only_derived_massless_cut_content() {
        let model = fermion_model();
        let options = GenerationOptions::default().with_graph_filter(
            GenerationFilter::PerturbativeOrders(BTreeMap::from([("QED".to_owned(), 1)])),
        );
        let unresolved = ResolvedProcess::unresolved_cut_content(&model, &options)
            .unwrap()
            .unwrap();
        assert_eq!(unresolved.maximum_multiplicity, 1);
        assert_eq!(unresolved.particles, BTreeSet::from([-1, 1, 22]));
        assert!(ResolvedProcess::matches_cut_content(
            &[1, 22],
            &[1],
            Some(&unresolved)
        ));
        assert!(!ResolvedProcess::matches_cut_content(
            &[1, 22, 22],
            &[1],
            Some(&unresolved)
        ));
        assert!(!ResolvedProcess::matches_cut_content(
            &[1, 25],
            &[1],
            Some(&unresolved)
        ));
    }

    #[test]
    fn exact_cut_rejects_species_only_false_positives() {
        let model = exact_cut_model();
        let interaction = |name: &str| ColoredNode::Interaction(name.to_owned());
        let external = |index, state, name, pdg| {
            ColoredNode::External(ExternalNode {
                index,
                state,
                particle: ParticleReference::new(name, pdg),
                symmetry_class: index + 1,
            })
        };
        let edge = |pdg| EdgeColor {
            pdg,
            direction: None,
        };

        let mut three_loop = Graph::new();
        let h0 = three_loop.add_node(interaction("V_78"));
        let h1 = three_loop.add_node(interaction("V_78"));
        let g0 = three_loop.add_node(interaction("V_76"));
        let p0 = three_loop.add_node(interaction("V_73"));
        let g1 = three_loop.add_node(interaction("V_76"));
        let p1 = three_loop.add_node(interaction("V_73"));
        let q0 = three_loop.add_node(interaction("V_98"));
        let q1 = three_loop.add_node(interaction("V_98"));
        let x0 = three_loop.add_node(external(0, ExternalState::Incoming, "e+", -11));
        let x1 = three_loop.add_node(external(1, ExternalState::Incoming, "e-", 11));
        let x2 = three_loop.add_node(external(2, ExternalState::Outgoing, "e+", -11));
        let x3 = three_loop.add_node(external(3, ExternalState::Outgoing, "e-", 11));
        three_loop.add_edge(h0, h1, false, edge(25)).unwrap();
        three_loop.add_edge(h0, p1, true, edge(5)).unwrap();
        three_loop.add_edge(h1, p0, true, edge(5)).unwrap();
        three_loop.add_edge(g0, g1, false, edge(21)).unwrap();
        three_loop.add_edge(g0, g1, true, edge(5)).unwrap();
        three_loop.add_edge(p0, g0, true, edge(5)).unwrap();
        three_loop.add_edge(p0, q1, false, edge(22)).unwrap();
        three_loop.add_edge(g1, h1, true, edge(5)).unwrap();
        three_loop.add_edge(p1, h0, true, edge(5)).unwrap();
        three_loop.add_edge(p1, q0, false, edge(22)).unwrap();
        three_loop.add_edge(q1, x0, true, edge(11)).unwrap();
        three_loop.add_edge(x1, q1, true, edge(11)).unwrap();
        three_loop.add_edge(q0, x2, true, edge(11)).unwrap();
        three_loop.add_edge(x3, q0, true, edge(11)).unwrap();

        let mut four_loop = Graph::new();
        let x0 = four_loop.add_node(external(0, ExternalState::Incoming, "e+", -11));
        let x1 = four_loop.add_node(external(1, ExternalState::Incoming, "e-", 11));
        let x2 = four_loop.add_node(external(2, ExternalState::Outgoing, "e+", -11));
        let x3 = four_loop.add_node(external(3, ExternalState::Outgoing, "e-", 11));
        let p0 = four_loop.add_node(interaction("V_73"));
        let p1 = four_loop.add_node(interaction("V_73"));
        let h0 = four_loop.add_node(interaction("V_78"));
        let h1 = four_loop.add_node(interaction("V_78"));
        let g0 = four_loop.add_node(interaction("V_76"));
        let g1 = four_loop.add_node(interaction("V_76"));
        let g2 = four_loop.add_node(interaction("V_76"));
        let g3 = four_loop.add_node(interaction("V_76"));
        let q0 = four_loop.add_node(interaction("V_98"));
        let q1 = four_loop.add_node(interaction("V_98"));
        four_loop.add_edge(q1, x0, true, edge(11)).unwrap();
        four_loop.add_edge(x1, q1, true, edge(11)).unwrap();
        four_loop.add_edge(p0, q1, false, edge(22)).unwrap();
        four_loop.add_edge(g3, p0, true, edge(5)).unwrap();
        four_loop.add_edge(p0, g2, true, edge(5)).unwrap();
        four_loop.add_edge(g2, g3, false, edge(21)).unwrap();
        four_loop.add_edge(h0, g3, true, edge(5)).unwrap();
        four_loop.add_edge(g2, h0, true, edge(5)).unwrap();
        four_loop.add_edge(h0, h1, false, edge(25)).unwrap();
        four_loop.add_edge(g0, h1, true, edge(5)).unwrap();
        four_loop.add_edge(h1, g1, true, edge(5)).unwrap();
        four_loop.add_edge(p1, g0, true, edge(5)).unwrap();
        four_loop.add_edge(g1, p1, true, edge(5)).unwrap();
        four_loop.add_edge(g0, g1, false, edge(21)).unwrap();
        four_loop.add_edge(p1, q0, false, edge(22)).unwrap();
        four_loop.add_edge(q0, x2, true, edge(11)).unwrap();
        four_loop.add_edge(x3, q0, true, edge(11)).unwrap();

        let resolved = ResolvedProcess::new(
            &model,
            &Process::cross_section([-11_i64, 11], [5_i64, -5, 25]),
        )
        .unwrap();
        // These graphs contain the requested species, but neither has an
        // s-channel cut with the exact final-state multiset after allowing only
        // the configured unresolved QCD content.
        for (graph, qcd_order) in [(three_loop, 1), (four_loop, 2)] {
            let options = GenerationOptions::default()
                .with_graph_filter(GenerationFilter::PerturbativeOrders(BTreeMap::from([(
                    "QCD".to_owned(),
                    qcd_order,
                )])))
                .with_graph_filter(GenerationFilter::BlobRange(1..=1))
                .with_graph_filter(GenerationFilter::SpectatorRange(0..=0));
            assert!(
                resolved
                    .cut_partitions(&graph, &model, &options)
                    .unwrap()
                    .is_empty()
            );
        }
    }

    #[test]
    fn sewn_filter_vetoes_factorization_created_by_pairwise_sewing() {
        let model = fermion_model();
        let process =
            ResolvedProcess::new(&model, &Process::cross_section([22_i64, 22], [22_i64])).unwrap();
        let mut graph = Graph::new();
        let externals: Vec<_> = (0..4)
            .map(|index| {
                graph.add_node(ColoredNode::External(ExternalNode {
                    index,
                    state: if index < 2 {
                        ExternalState::Incoming
                    } else {
                        ExternalState::Outgoing
                    },
                    particle: ParticleReference::new("a", 22),
                    symmetry_class: index + 1,
                }))
            })
            .collect();
        let center = graph.add_node(ColoredNode::Interaction("V".to_owned()));
        let photon = EdgeColor {
            pdg: 22,
            direction: None,
        };
        for external in &externals {
            graph.add_edge(*external, center, false, photon).unwrap();
        }

        let options = GenerationOptions::default()
            .with_graph_filter(GenerationFilter::Sewn(SewnFilterOptions::default()));
        assert!(process.passes_sewn_filter(&graph, &options));

        let mut factorized = Graph::new();
        let externals: Vec<_> = (0..4)
            .map(|index| {
                factorized.add_node(ColoredNode::External(ExternalNode {
                    index,
                    state: if index < 2 {
                        ExternalState::Incoming
                    } else {
                        ExternalState::Outgoing
                    },
                    particle: ParticleReference::new("a", 22),
                    symmetry_class: index + 1,
                }))
            })
            .collect();
        let left = factorized.add_node(ColoredNode::Interaction("V".to_owned()));
        let right = factorized.add_node(ColoredNode::Interaction("V".to_owned()));
        for external in [externals[0], externals[2]] {
            factorized.add_edge(external, left, false, photon).unwrap();
        }
        for external in [externals[1], externals[3]] {
            factorized.add_edge(external, right, false, photon).unwrap();
        }
        factorized.add_edge(left, right, false, photon).unwrap();
        assert!(!process.passes_sewn_filter(&factorized, &options));
        assert!(process.passes_sewn_filter(
            &factorized,
            &GenerationOptions::default().with_graph_filter(GenerationFilter::Sewn(
                SewnFilterOptions {
                    filter_tadpoles: false,
                }
            ))
        ));
    }

    #[test]
    fn amplitude_fermions_require_explicit_symmetrization_opt_in() {
        let model = fermion_model();
        let classes = |process: Process| {
            ResolvedProcess::new(&model, &process)
                .unwrap()
                .external_edges(&model)
                .unwrap()
                .into_iter()
                .map(|(node, _)| node.external.unwrap().symmetry_class)
                .collect::<Vec<_>>()
        };
        let process =
            Process::amplitude([1_i64, -1, 22], Vec::<i64>::new()).symmetrize_initial(true);
        assert_eq!(classes(process.clone()), vec![1, 2, 0]);
        assert_eq!(
            classes(process.symmetrize_external_fermions(true)),
            vec![0, 0, 0]
        );

        let cross_section = Process::cross_section([1_i64, -1], [22_i64]).symmetrize_initial(true);
        assert_eq!(classes(cross_section), vec![0, 0, 3, 4]);
    }

    #[test]
    fn numerator_construction_localizes_internal_ufo_data() {
        let mut definition = fermion_model().into_definition();
        definition
            .particles
            .iter_mut()
            .find(|particle| particle.name == "f")
            .unwrap()
            .propagator = Some("custom_f_prop".to_owned());
        definition.propagators.push(Propagator {
            name: "custom_f_prop".to_owned(),
            particle: "f".to_owned(),
            numerator: "UFO::{}::P(UFO::{}::idx(1,1))*UFO::{}::P(UFO::{}::idx(1,2))*UFO::{}::Metric(UFO::{}::idx(1,1),UFO::{}::dummy(1))*UFO::{}::Metric(UFO::{}::dummy(1),UFO::{}::idx(1,2))".to_owned(),
            denominator: "P^2".to_owned(),
        });
        let generator = Generator::new(Model::new(definition).unwrap());

        let mut graph = Graph::new();
        let interaction = graph.add_node(ColoredNode::Interaction("V".to_owned()));
        let external = graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Outgoing,
            particle: ParticleReference::new("a", 22),
            symmetry_class: 1,
        }));
        graph
            .add_edge(
                interaction,
                interaction,
                true,
                EdgeColor {
                    pdg: 1,
                    direction: None,
                },
            )
            .unwrap();
        graph
            .add_edge(
                interaction,
                external,
                false,
                EdgeColor {
                    pdg: 22,
                    direction: None,
                },
            )
            .unwrap();

        let diagram = generator
            .to_diagram(
                "localized".to_owned(),
                ColoredTopology {
                    graph,
                    symmetry: 1,
                    multiplicity: 1,
                    cut_partitions: Vec::new(),
                },
                1,
                FermionSigns::default(),
            )
            .unwrap();
        let internal_numerator = diagram
            .edges()
            .find(|(id, _, _)| id.0 == 0)
            .unwrap()
            .2
            .numerator
            .as_deref()
            .unwrap();
        assert!(internal_numerator.contains("FeynKit::Momentum"));
        assert!(internal_numerator.contains("FeynKit::SourceIndex"));
        assert!(internal_numerator.contains("FeynKit::SinkIndex"));
        assert!(internal_numerator.contains("FeynKit::EdgeDummy"));
        assert!(!internal_numerator.contains("idx("));
        assert!(!internal_numerator.contains("Slash(P)"));

        let external_numerator = diagram
            .edges()
            .find(|(id, _, _)| id.0 == 1)
            .unwrap()
            .2
            .numerator
            .as_deref();
        assert_eq!(external_numerator, Some("1"));

        let vertex_numerator = diagram
            .vertex(feynkit_graph::VertexId(interaction))
            .unwrap()
            .numerator
            .as_deref()
            .unwrap();
        assert!(vertex_numerator.contains("FeynKit::SourceIndex"));
        assert!(vertex_numerator.contains("FeynKit::SinkIndex"));
    }
}
