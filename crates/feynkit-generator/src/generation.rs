// The case table of the generation oracle tests, shared by the unit tests that run its workloads.
#[cfg(test)]
#[path = "../tests/generation_oracle/cases.rs"]
#[allow(dead_code)]
mod cases;
mod enumeration;
mod graph;
#[cfg(test)]
mod projector_equivalence;
#[cfg(test)]
mod search_census;
mod special_topologies;

use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fmt,
    hash::Hash,
    sync::{
        Arc, Mutex, RwLock,
        atomic::{AtomicUsize, Ordering},
    },
};

use feynkit_graph::{
    DiagramCut, DiagramCutSide, DiagramEdge, DiagramEndpoint, DiagramError, DiagramHalfEdge,
    DiagramId, DiagramThresholdCandidate, DiagramVertex, DotDiagram, EdgeId, ExternalLeg,
    ExternalState, FeynmanDiagram, FeynmanDiagramBuilder, VertexId, VertexSlot,
};
use feynkit_model::{Model, ModelError, Particle, ParticleId, VertexRule, VertexRuleId};
use idenso::{color::CS, dirac::AGS, epsilon::EPSILON_SYMBOL, representations::Bispinor};
use linnet::half_edge::{
    HedgeGraph, NodeIndex,
    builder::HedgeGraphBuilder,
    canon::ignore,
    generate::GraphGenerator,
    involution::{Flow, Orientation},
    subgraph::{InternalSubGraph, SuBitGraph, SubSetLike},
};
use rayon::ThreadPoolBuildError;
#[cfg(not(target_arch = "wasm32"))]
use rayon::{ThreadPoolBuilder, prelude::*};
use serde::{Deserialize, Serialize};
use spenso::{
    network::library::symbolic::ETS,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    symbol,
};
use thiserror::Error;

use feynkit_amplitude::ColorRepresentation;

use crate::{
    FilterScope, GenerationFilter, GenerationFilterKind, GenerationOptions, GenerationOutput,
    GenerationType, GroupingError, NumeratorGrouping, ParticleSelector, Process, ProcessError,
    SelectorError, SewnFilterOptions, grouping, momentum_symbol,
};
pub use enumeration::{ExternalLabel, GenerationView};
use enumeration::{
    ParticleStub, ProcessExternals, Pruner, VertexSignatures, closed_loops, incoming_stub,
    resolved_vertex,
};
pub use graph::ExternalNode;
use graph::{Endpoint, GeneratedEdge, GeneratedGraph, PlainNodes, PlainVertex};

#[derive(Debug, Error)]
pub enum GenerationError {
    #[error("name selection and named momentum bases require named generation output")]
    NamesRequireNamedOutput,
    #[error(transparent)]
    Process(#[from] ProcessError),
    #[error(transparent)]
    Amplitude(#[from] Box<feynkit_amplitude::AmplitudeError>),
    #[error("generation was cancelled; cannot construct a complete amplitude")]
    IncompleteAmplitude,
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error(transparent)]
    Selector(#[from] SelectorError),
    #[error(transparent)]
    Diagram(#[from] DiagramError),
    #[error("failed to build generation thread pool: {0}")]
    ThreadPool(#[from] ThreadPoolBuildError),
    #[error("graph automorphism factor cannot be represented by u64")]
    SymmetryFactor,
    #[error(transparent)]
    Enumeration(#[from] linnet::half_edge::generate::GenerationError),
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
    #[error(
        "symmetrize_final does not apply to cross sections: identical final-state particles are cut edges, already counted by the graph automorphisms; remove the option"
    )]
    CrossSectionFinalStateSymmetrization,
    #[error(
        "symmetrize_initial exchanges the two identical beams of a 2 -> X cross section, but this process has {incoming} incoming particles with a repeated species, which no Lorentz transformation permutes"
    )]
    UnsupportedBeamSymmetrization { incoming: usize },
    #[error("generation thread count must be greater than zero, received {0}")]
    InvalidThreadCount(usize),
    #[error("numerator grouping requires at least one numerical sample")]
    InvalidNumericalSampleCount,
    #[error("arithmetic overflow while computing {0}")]
    ArithmeticOverflow(&'static str),
    #[error("forced cut edge set {edges:?} does not match any generated physical cut")]
    UnknownForcedCut { edges: Vec<usize> },
    #[error(transparent)]
    Grouping(#[from] GroupingError),
    #[error("invalid generation grouping provenance: {0}")]
    InvalidGroups(String),
    #[error(
        "fermion flow branches at vertex {vertex} with degree {degree}; branching fermion chains are not supported"
    )]
    UnsupportedFermionBranching { vertex: usize, degree: usize },
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
    #[error("the {owner} numerator has an invalid UFO tensor {tensor}: {message}")]
    InvalidNumeratorTensor {
        owner: String,
        tensor: String,
        message: String,
    },
    #[error("the {owner} numerator uses unsupported UFO tensor {tensor}")]
    UnsupportedNumeratorTensor { owner: String, tensor: String },
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
    #[error(
        "could not identify a unique CP-conjugate rule for interaction '{interaction}' on vertex {vertex}; compatible rules: {candidates:?}"
    )]
    InvalidCpInteractionRule {
        vertex: usize,
        interaction: String,
        candidates: Vec<String>,
    },
}

/// Source metadata requested by named generation output.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceIdentity {
    pub id: DiagramId,
    pub name: String,
}

#[derive(Debug, Clone)]
pub struct GroupMember {
    /// Stable generated-order index assigned before zero-numerator removal.
    pub source_diagram: usize,
    /// Requested source identity; absent for minimal output.
    pub source_identity: Option<SourceIdentity>,
    /// Index into [`GenerationResult::diagrams`] after zero-numerator removal.
    pub diagram: usize,
    /// The numerator of this member divided by the numerator of the master.
    pub ratio: Atom,
    /// Numerator-independent factor contributed by this source diagram.
    pub overall_factor: Atom,
}

#[derive(Debug, Clone)]
pub struct DiagramGroup {
    /// Post-zero-removal index of the numerator reference diagram.
    pub master: usize,
    /// Deterministically ordered members, including the master with ratio one.
    pub members: Vec<GroupMember>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenerationReport {
    /// Complete untyped topologies accepted during enumeration. Loop orders that the coupling
    /// orders rule out (no multiset of vertex rules has the degrees and orders they need), and
    /// requests with external charges or parities that no vertex rule balances or with an empty
    /// coupling-order or fermion-loop range, are not enumerated and add none.
    pub topology_count: usize,
    /// Enumerated diagrams: topologies with fields and vertex rules, before the cut and
    /// charge-conjugation filters that need their canonical frame.
    pub interaction_assignment_count: usize,
    /// Number of diagrams retained after exact zero-numerator removal.
    pub retained_count: usize,
    pub zero_numerator_count: usize,
    pub completed: bool,
}

#[derive(Debug, Clone)]
pub struct GenerationResult {
    /// Final retained master diagrams in deterministic order.
    pub diagrams: Vec<FeynmanDiagram>,
    pub groups: Vec<DiagramGroup>,
    pub report: GenerationReport,
}

impl GenerationResult {
    /// Convert complete unsewn diagrams into a coherent symbolic amplitude.
    pub fn into_amplitude(
        self,
        options: feynkit_amplitude::AmplitudeOptions,
    ) -> Result<feynkit_amplitude::Amplitude, GenerationError> {
        if !self.report.completed {
            return Err(GenerationError::IncompleteAmplitude);
        }
        Ok(
            feynkit_amplitude::Amplitude::new(self.diagrams.into_iter().map(Arc::new), options)
                .map_err(Box::new)?,
        )
    }

    /// Validate collapsed group provenance against the finalized master list.
    pub fn validate_groups(&self) -> Result<(), GenerationError> {
        if self.groups.len() != self.diagrams.len() {
            return Err(GenerationError::InvalidGroups(format!(
                "{} groups do not cover {} collapsed diagrams",
                self.groups.len(),
                self.diagrams.len()
            )));
        }
        let mut masters = BTreeSet::new();
        for (group_index, group) in self.groups.iter().enumerate() {
            if group.master >= self.diagrams.len() || !masters.insert(group.master) {
                return Err(GenerationError::InvalidGroups(format!(
                    "group {group_index} has invalid or duplicate master {}",
                    group.master
                )));
            }
            if group.members.is_empty()
                || group
                    .members
                    .iter()
                    .any(|member| member.diagram != group.master)
            {
                return Err(GenerationError::InvalidGroups(format!(
                    "group {group_index} has no members or a member assigned to another master"
                )));
            }
            self.diagrams[group.master].validate()?;
        }
        if masters != (0..self.diagrams.len()).collect::<BTreeSet<_>>() {
            return Err(GenerationError::InvalidGroups(
                "collapsed diagram masters are not covered exactly once".to_owned(),
            ));
        }
        Ok(())
    }
}

/// Extra massless particles allowed by perturbative unresolved-cut selection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnresolvedCutContent {
    pub particles: BTreeSet<ParticleId>,
    pub maximum_multiplicity: usize,
}

/// Resolve the canonical unresolved-cut policy stored in generation options.
pub fn unresolved_cut_content(
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
        let orders = rule.coupling_orders(model);
        if !required.keys().any(|name| orders.contains_key(name)) {
            continue;
        }
        for particle in &rule.particles {
            if model.particle_is_massless(*particle) {
                particles.insert(*particle);
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct NumeratorHalfEdge {
    edge: usize,
    // Native half-edge IDs need not equal 2 * edge + endpoint after sewing.
    hedge: usize,
    flow: Flow,
    spin: i64,
    color: i64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum NumeratorOwner {
    Vertex(usize),
    Edge(usize),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NumeratorSector {
    Spin,
    Color,
}

impl fmt::Display for NumeratorOwner {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Vertex(vertex) => write!(formatter, "vertex {vertex}"),
            Self::Edge(edge) => write!(formatter, "edge {edge}"),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum NumeratorFrame {
    Diagram,
    Interaction,
}

struct NumeratorInstantiation<'a> {
    owner: NumeratorOwner,
    legs: &'a [NumeratorHalfEdge],
    frame: NumeratorFrame,
}

impl<'legs> NumeratorInstantiation<'legs> {
    fn new(owner: NumeratorOwner, legs: &'legs [NumeratorHalfEdge]) -> Self {
        Self {
            owner,
            legs,
            frame: NumeratorFrame::Diagram,
        }
    }

    fn dummy_symbol(&self) -> symbolica::atom::Symbol {
        if self.frame == NumeratorFrame::Interaction {
            symbol!("feynkit_generator_template::dummy")
        } else {
            match self.owner {
                NumeratorOwner::Vertex(_) => feynkit_graph::symbols::vertex_index(),
                NumeratorOwner::Edge(_) => feynkit_graph::symbols::edge_index(),
            }
        }
    }

    fn momentum_symbol(&self) -> symbolica::atom::Symbol {
        if self.frame == NumeratorFrame::Interaction {
            symbol!("feynkit_generator_template::momentum")
        } else {
            momentum_symbol()
        }
    }

    // UFO `idx(shift, leg)` and `dummy(index)` values are local to one
    // propagator or interaction. Localize them before exposing the annotation
    // so equal endpoint indices contract and unrelated dummy indices cannot.
    fn instantiate(
        &self,
        expression: &Atom,
        sector: NumeratorSector,
    ) -> Result<Atom, GenerationError> {
        let mut atom = expression.clone();
        let mut error = None;
        atom = atom.replace_map(|term, _, out| {
            if error.is_some() {
                return;
            }
            match self.localize_momentum(term) {
                Ok(Some(replacement)) => {
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
                    **out = replacement;
                }
                Ok(None) => {}
                Err(localization_error) => error = Some(localization_error),
            }
        });
        if let Some(error) = error {
            return Err(error);
        }

        atom = match sector {
            NumeratorSector::Spin => self.lower_spin(atom)?,
            NumeratorSector::Color => self.lower_color(atom)?,
        };

        Ok(atom)
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
                        if self.frame == NumeratorFrame::Interaction {
                            symbol!("feynkit_generator_template::orientation").call((leg.edge,))
                                * momentum
                        } else if leg.flow == Flow::Source {
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

    fn lower_spin(&self, mut atom: Atom) -> Result<Atom, GenerationError> {
        // Force the canonical Idenso/Spenso registrations before constructing
        // any of their heads. Symbolica attributes cannot be added after a
        // bare symbol with the same name has already been interned.
        let _ = (
            AGS.gamma,
            AGS.gamma5,
            AGS.charge_conjugation,
            AGS.projm,
            AGS.projp,
            AGS.sigma,
            ETS.metric,
            *EPSILON_SYMBOL,
        );

        let mut next_dummy = self.maximum_local_dummy(&atom)?;
        let mut error = None;
        atom = atom.replace_map(|term, _, out| {
            if error.is_some() {
                return;
            }
            let AtomView::Fun(function) = term else {
                return;
            };
            let symbol = function.get_symbol();
            let name = symbol.get_stripped_name();
            if !Self::is_model_symbol(symbol, name) {
                return;
            }
            let arguments = function
                .iter()
                .map(|argument| argument.to_owned())
                .collect::<Vec<_>>();
            let lowered = match name {
                "Identity" => self.exact_arguments(term, &arguments, 2).map(|args| {
                    ETS.metric(
                        self.bispinor(args[0].clone()),
                        self.bispinor(args[1].clone()),
                    )
                }),
                "IdentityL" | "Metric" => self.exact_arguments(term, &arguments, 2).map(|args| {
                    ETS.metric(
                        self.minkowski(args[0].clone()),
                        self.minkowski(args[1].clone()),
                    )
                }),
                "Gamma" => self.exact_arguments(term, &arguments, 3).map(|args| {
                    FunctionBuilder::new(AGS.gamma)
                        .add_arg(self.bispinor(args[1].clone()))
                        .add_arg(self.bispinor(args[2].clone()))
                        .add_arg(self.minkowski(args[0].clone()))
                        .finish()
                }),
                "Gamma5" => self.exact_arguments(term, &arguments, 2).map(|args| {
                    FunctionBuilder::new(AGS.gamma5)
                        .add_arg(self.bispinor(args[0].clone()))
                        .add_arg(self.bispinor(args[1].clone()))
                        .finish()
                }),
                "ProjM" => self.exact_arguments(term, &arguments, 2).map(|args| {
                    FunctionBuilder::new(AGS.projm)
                        .add_arg(self.bispinor(args[0].clone()))
                        .add_arg(self.bispinor(args[1].clone()))
                        .finish()
                }),
                "ProjP" => self.exact_arguments(term, &arguments, 2).map(|args| {
                    FunctionBuilder::new(AGS.projp)
                        .add_arg(self.bispinor(args[0].clone()))
                        .add_arg(self.bispinor(args[1].clone()))
                        .finish()
                }),
                "Sigma" => self.exact_arguments(term, &arguments, 4).map(|args| {
                    FunctionBuilder::new(AGS.sigma)
                        .add_arg(self.bispinor(args[2].clone()))
                        .add_arg(self.bispinor(args[3].clone()))
                        .add_arg(self.minkowski(args[0].clone()))
                        .add_arg(self.minkowski(args[1].clone()))
                        .finish()
                }),
                "PSlash" => {
                    self.exact_arguments(term, &arguments, 3).and_then(|args| {
                        next_dummy = next_dummy.checked_add(1).ok_or(
                            GenerationError::ArithmeticOverflow("a symbolic slash dummy index"),
                        )?;
                        let index = self.minkowski(self.dummy(next_dummy)?);
                        let momentum = self.index_momentum(args[2].as_view(), index.clone())?;
                        Ok(FunctionBuilder::new(AGS.gamma)
                            .add_arg(self.bispinor(args[0].clone()))
                            .add_arg(self.bispinor(args[1].clone()))
                            .add_arg(index)
                            .finish()
                            * momentum)
                    })
                }
                "Epsilon" => {
                    if arguments.is_empty() {
                        Err(self.invalid_tensor(term, "expected at least one Lorentz index"))
                    } else {
                        Ok(arguments
                            .into_iter()
                            .fold(
                                FunctionBuilder::new(*EPSILON_SYMBOL),
                                |builder, argument| builder.add_arg(self.minkowski(argument)),
                            )
                            .finish())
                    }
                }
                "C" => self.exact_arguments(term, &arguments, 2).map(|args| {
                    FunctionBuilder::new(AGS.charge_conjugation)
                        .add_arg(self.bispinor(args[0].clone()))
                        .add_arg(self.bispinor(args[1].clone()))
                        .finish()
                }),
                "T" | "f" | "d" | "EpsilonBar" | "T6" | "K6" | "K6Bar" => {
                    Err(self.invalid_tensor(term, "a color tensor appeared in a Lorentz structure"))
                }
                _ => return,
            };
            match lowered {
                Ok(replacement) => **out = replacement,
                Err(lowering_error) => error = Some(lowering_error),
            }
        });
        if let Some(error) = error {
            return Err(error);
        }

        // P(mu, leg) has already become gammalooprs::Q(edge, mu).
        // Attach the Minkowski representation after all UFO indices have been
        // localized. PSlash momenta already carry the freshly allocated slot.
        atom = atom.replace_map(|term, _, out| {
            let AtomView::Fun(function) = term else {
                return;
            };
            if function.get_symbol() != self.momentum_symbol() || function.get_nargs() != 2 {
                return;
            }
            let mut arguments = function.iter();
            let edge = arguments.next().expect("the arity was checked");
            let index = arguments.next().expect("the arity was checked");
            if Self::is_representation(index, "mink") {
                return;
            }
            **out = FunctionBuilder::new(self.momentum_symbol())
                .add_arg(edge)
                .add_arg(self.minkowski(index.to_owned()))
                .finish();
        });
        self.reject_residual_tensors(&atom)?;
        Ok(atom)
    }

    fn lower_color(&self, mut atom: Atom) -> Result<Atom, GenerationError> {
        let _ = (CS.t, CS.f, ETS.metric, *EPSILON_SYMBOL);
        atom = atom.expand();
        if let AtomView::Add(terms) = atom.as_view() {
            return terms.iter().try_fold(Atom::Zero, |sum, term| {
                Ok(sum + self.lower_color(term.to_owned())?)
            });
        }
        let dummy_symbol = self.dummy_symbol();
        // Contract identities touching a local dummy before assigning oriented
        // slots. Each step removes one identity, including chains of identities;
        // the remaining tensors or external legs determine its representation.
        loop {
            let factors = match atom.as_view() {
                AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                factor => vec![factor],
            };
            let mut contraction = None;
            for term in factors {
                let AtomView::Fun(function) = term else {
                    continue;
                };
                if !Self::is_model_symbol(function.get_symbol(), "Identity")
                    || function.get_nargs() != 2
                {
                    continue;
                }
                let args = function.iter().collect::<Vec<_>>();
                if args[0] == args[1] {
                    continue;
                }
                for (left, right) in [(args[0], args[1]), (args[1], args[0])] {
                    if !matches!(left, AtomView::Fun(index) if index.get_symbol() == dummy_symbol) {
                        continue;
                    }
                    let mut occurrences = 0;
                    let mut powered = false;
                    atom.visitor(&mut |index| {
                        if let AtomView::Pow(_) = index {
                            index.visitor(&mut |nested| {
                                powered |= nested == left;
                                true
                            });
                            return false;
                        }
                        occurrences += usize::from(index == left);
                        true
                    });
                    if occurrences == 2 && !powered {
                        contraction = Some((term.to_owned(), left.to_owned(), right.to_owned()));
                        break;
                    }
                }
                if contraction.is_some() {
                    break;
                }
            }
            let Some((identity, dummy, other)) = contraction else {
                break;
            };
            atom = (atom / identity).replace(dummy.to_pattern()).with(other);
        }
        let mut representations = HashMap::<String, ColorRepresentation>::new();
        for leg in self.legs {
            if let Some(representation) = ColorRepresentation::from_ufo(leg.color) {
                let index = self.index(*leg, 1)?;
                self.assign_color_representation(
                    &mut representations,
                    index.as_view(),
                    representation,
                )?;
            }
        }

        let mut identities = Vec::<(Atom, Atom)>::new();
        let mut error = None;
        atom = atom.replace_map(|term, _, _| {
            if error.is_some() {
                return;
            }
            let AtomView::Fun(function) = term else {
                return;
            };
            let symbol = function.get_symbol();
            let name = symbol.get_stripped_name();
            if !Self::is_model_symbol(symbol, name) {
                return;
            }
            let arguments = function
                .iter()
                .map(|argument| argument.to_owned())
                .collect::<Vec<_>>();
            let result = match name {
                "Identity" => self.exact_arguments(term, &arguments, 2).map(|args| {
                    identities.push((args[0].clone(), args[1].clone()));
                }),
                "T" => self.exact_arguments(term, &arguments, 3).and_then(|args| {
                    self.assign_color_representation(
                        &mut representations,
                        args[0].as_view(),
                        ColorRepresentation::Adjoint,
                    )?;
                    self.assign_color_representation(
                        &mut representations,
                        args[1].as_view(),
                        ColorRepresentation::Fundamental,
                    )?;
                    self.assign_color_representation(
                        &mut representations,
                        args[2].as_view(),
                        ColorRepresentation::AntiFundamental,
                    )
                }),
                "f" | "d" => self.exact_arguments(term, &arguments, 3).and_then(|args| {
                    for argument in args {
                        self.assign_color_representation(
                            &mut representations,
                            argument.as_view(),
                            ColorRepresentation::Adjoint,
                        )?;
                    }
                    Ok(())
                }),
                "Epsilon" | "EpsilonBar" => {
                    self.exact_arguments(term, &arguments, 3).and_then(|args| {
                        let representation = if name == "Epsilon" {
                            ColorRepresentation::Fundamental
                        } else {
                            ColorRepresentation::AntiFundamental
                        };
                        for argument in args {
                            self.assign_color_representation(
                                &mut representations,
                                argument.as_view(),
                                representation,
                            )?;
                        }
                        Ok(())
                    })
                }
                "T6" | "K6" | "K6Bar" => Err(self.unsupported_tensor(term)),
                "IdentityL" | "Gamma" | "Gamma5" | "ProjM" | "ProjP" | "Sigma" | "C" | "Metric"
                | "PSlash" => {
                    Err(self.invalid_tensor(term, "a Lorentz tensor appeared in a color structure"))
                }
                _ => return,
            };
            if let Err(lowering_error) = result {
                error = Some(lowering_error);
            }
        });
        if let Some(error) = error {
            return Err(error);
        }

        loop {
            let mut changed = false;
            for (left, right) in &identities {
                let left_representation = representations
                    .get(&Self::index_key(left.as_view()))
                    .copied();
                let right_representation = representations
                    .get(&Self::index_key(right.as_view()))
                    .copied();
                match (left_representation, right_representation) {
                    (Some(left_rep), Some(right_rep)) if left_rep.dual() != right_rep => {
                        return Err(self.invalid_tensor(
                            FunctionBuilder::new(symbol!("UFO::Identity"))
                                .add_arg(left)
                                .add_arg(right)
                                .finish()
                                .as_view(),
                            "the two color identity indices do not carry dual representations",
                        ));
                    }
                    (Some(_), Some(_)) => {}
                    (Some(left), None) => {
                        self.assign_color_representation(
                            &mut representations,
                            right.as_view(),
                            left.dual(),
                        )?;
                        changed = true;
                    }
                    (None, Some(right)) => {
                        self.assign_color_representation(
                            &mut representations,
                            left.as_view(),
                            right.dual(),
                        )?;
                        changed = true;
                    }
                    (None, None) => {}
                }
            }
            if !changed {
                break;
            }
        }
        for (left, right) in &identities {
            if !representations.contains_key(&Self::index_key(left.as_view()))
                || !representations.contains_key(&Self::index_key(right.as_view()))
            {
                return Err(self.invalid_tensor(
                    FunctionBuilder::new(symbol!("UFO::Identity"))
                        .add_arg(left)
                        .add_arg(right)
                        .finish()
                        .as_view(),
                    "its color representation cannot be inferred",
                ));
            }
        }

        error = None;
        atom = atom.replace_map(|term, _, out| {
            if error.is_some() {
                return;
            }
            let AtomView::Fun(function) = term else {
                return;
            };
            let symbol = function.get_symbol();
            let name = symbol.get_stripped_name();
            if !Self::is_model_symbol(symbol, name) {
                return;
            }
            let arguments = function
                .iter()
                .map(|argument| argument.to_owned())
                .collect::<Vec<_>>();
            let lowered = match name {
                "Identity" => self.exact_arguments(term, &arguments, 2).map(|args| {
                    let left = representations[&Self::index_key(args[0].as_view())];
                    let right = representations[&Self::index_key(args[1].as_view())];
                    ETS.metric(
                        left.representation().to_symbolic([args[0].clone()]),
                        right.representation().to_symbolic([args[1].clone()]),
                    )
                }),
                "T" => self.exact_arguments(term, &arguments, 3).map(|args| {
                    FunctionBuilder::new(CS.t)
                        .add_arg(
                            ColorRepresentation::Adjoint
                                .representation()
                                .to_symbolic([args[0].clone()]),
                        )
                        .add_arg(
                            ColorRepresentation::Fundamental
                                .representation()
                                .to_symbolic([args[1].clone()]),
                        )
                        .add_arg(
                            ColorRepresentation::AntiFundamental
                                .representation()
                                .to_symbolic([args[2].clone()]),
                        )
                        .finish()
                }),
                "f" | "d" => self.exact_arguments(term, &arguments, 3).map(|args| {
                    let adjoint_indices = args.iter().map(|argument| {
                        ColorRepresentation::Adjoint
                            .representation()
                            .to_symbolic([argument.clone()])
                    });
                    if name == "d" {
                        // UFO d^{abc} = 2 Tr({T^a,T^b} T^c), whereas Idenso
                        // stores the normalized symmetric generator trace.
                        Atom::num(4)
                            * CS.symmetric_generator_trace(
                                ColorRepresentation::Fundamental
                                    .representation()
                                    .to_symbolic([]),
                                adjoint_indices,
                            )
                    } else {
                        adjoint_indices
                            .fold(FunctionBuilder::new(CS.f), |builder, argument| {
                                builder.add_arg(argument)
                            })
                            .finish()
                    }
                }),
                "Epsilon" | "EpsilonBar" => self.exact_arguments(term, &arguments, 3).map(|args| {
                    let representation = if name == "Epsilon" {
                        ColorRepresentation::Fundamental
                    } else {
                        ColorRepresentation::AntiFundamental
                    };
                    idenso::epsilon::epsilon(args.iter().map(|argument| {
                        representation
                            .representation()
                            .to_symbolic([argument.clone()])
                    }))
                }),
                _ => return,
            };
            match lowered {
                Ok(replacement) => **out = replacement,
                Err(lowering_error) => error = Some(lowering_error),
            }
        });
        if let Some(error) = error {
            return Err(error);
        }
        self.reject_residual_tensors(&atom)?;
        Ok(atom)
    }

    fn edge_color_identity(&self, color: i64) -> Result<Option<Atom>, GenerationError> {
        let Some(source_representation) = ColorRepresentation::from_ufo(color) else {
            return Ok(None);
        };
        let source = self
            .legs
            .iter()
            .copied()
            .find(|leg| leg.flow == Flow::Source)
            .ok_or_else(|| GenerationError::MissingNumeratorMomentum {
                owner: self.owner.to_string(),
            })?;
        let sink = self
            .legs
            .iter()
            .copied()
            .find(|leg| leg.flow == Flow::Sink)
            .ok_or_else(|| GenerationError::MissingNumeratorMomentum {
                owner: self.owner.to_string(),
            })?;
        let identity = ETS.metric(
            source_representation
                .representation()
                .to_symbolic([self.index(source, 1)?]),
            source_representation
                .dual()
                .representation()
                .to_symbolic([self.index(sink, 1)?]),
        );
        Ok(Some(identity))
    }

    fn maximum_local_dummy(&self, atom: &Atom) -> Result<usize, GenerationError> {
        let mut maximum = 0;
        let mut error = None;
        let _ = atom.replace_map(|term, _, _| {
            if error.is_some() {
                return;
            }
            let AtomView::Fun(function) = term else {
                return;
            };
            let expected = self.dummy_symbol();
            if function.get_symbol() != expected || function.get_nargs() != 2 {
                return;
            }
            let Some(local) = function.iter().nth(1) else {
                return;
            };
            match usize::try_from(local) {
                Ok(local) => maximum = maximum.max(local),
                Err(_) => error = Some(self.invalid_index(term)),
            }
        });
        error.map_or(Ok(maximum), Err)
    }

    fn exact_arguments<'a>(
        &self,
        tensor: AtomView<'_>,
        arguments: &'a [Atom],
        expected: usize,
    ) -> Result<&'a [Atom], GenerationError> {
        if arguments.len() == expected {
            Ok(arguments)
        } else {
            Err(self.invalid_tensor(
                tensor,
                &format!("expected {expected} arguments, found {}", arguments.len()),
            ))
        }
    }

    fn index_momentum(&self, momentum: AtomView<'_>, index: Atom) -> Result<Atom, GenerationError> {
        let AtomView::Fun(function) = momentum else {
            return Err(self.invalid_tensor(momentum, "expected a localized momentum"));
        };
        if function.get_symbol() != self.momentum_symbol() || function.get_nargs() != 1 {
            return Err(self.invalid_tensor(momentum, "expected a rank-one localized momentum"));
        }
        Ok(FunctionBuilder::new(self.momentum_symbol())
            .add_arg(function.iter().next().expect("the arity was checked"))
            .add_arg(index)
            .finish())
    }

    fn bispinor(&self, index: Atom) -> Atom {
        Bispinor {}.new_rep(4).to_symbolic([index])
    }

    fn minkowski(&self, index: Atom) -> Atom {
        Minkowski {}.new_rep(4).to_symbolic([index])
    }

    fn is_representation(index: AtomView<'_>, name: &str) -> bool {
        matches!(
            index,
            AtomView::Fun(function)
                if function.get_symbol().get_namespace() == "spenso"
                    && function.get_symbol().get_stripped_name() == name
                    && matches!(function.get_nargs(), 1 | 2)
        )
    }

    fn assign_color_representation(
        &self,
        representations: &mut HashMap<String, ColorRepresentation>,
        index: AtomView<'_>,
        representation: ColorRepresentation,
    ) -> Result<(), GenerationError> {
        let key = Self::index_key(index);
        let local_dummy = matches!(index, AtomView::Fun(function)
            if [feynkit_graph::symbols::vertex_index(), feynkit_graph::symbols::edge_index(), self.dummy_symbol()]
                .contains(&function.get_symbol()));
        match representations.insert(key, representation) {
            // A contracted fundamental index joins dual slots. External ports
            // retain the exact orientation fixed by their particle record.
            Some(previous) if previous != representation
                && !(local_dummy && previous == representation.dual()) => Err(self.invalid_tensor(
                index,
                &format!(
                    "index is used as both {previous:?} and {representation:?} color representations"
                ),
            )),
            _ => Ok(()),
        }
    }

    fn index_key(index: AtomView<'_>) -> String {
        index.to_canonical_string()
    }

    fn reject_residual_tensors(&self, atom: &Atom) -> Result<(), GenerationError> {
        let mut residual = None;
        let _ = atom.replace_map(|term, _, _| {
            if residual.is_some() {
                return;
            }
            let AtomView::Fun(function) = term else {
                return;
            };
            let symbol = function.get_symbol();
            let name = symbol.get_stripped_name();
            if Self::is_model_symbol(symbol, name)
                && matches!(
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
                        | "PSlash"
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
                residual = Some(self.unsupported_tensor(term));
            }
        });
        residual.map_or(Ok(()), Err)
    }

    fn invalid_tensor(&self, tensor: AtomView<'_>, message: &str) -> GenerationError {
        GenerationError::InvalidNumeratorTensor {
            owner: self.owner.to_string(),
            tensor: tensor.to_string(),
            message: message.to_owned(),
        }
    }

    fn unsupported_tensor(&self, tensor: AtomView<'_>) -> GenerationError {
        GenerationError::UnsupportedNumeratorTensor {
            owner: self.owner.to_string(),
            tensor: tensor.to_string(),
        }
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
        let hedge = i64::try_from(leg.hedge)
            .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic half-edge identifier"))?;
        let symbol = if self.frame == NumeratorFrame::Interaction {
            symbol!("feynkit_generator_template::hedge")
        } else {
            feynkit_graph::symbols::hedge_index()
        };
        Ok(symbol.call((hedge, shift)))
    }

    fn dummy(&self, local: usize) -> Result<Atom, GenerationError> {
        let (symbol, owner) = match self.owner {
            NumeratorOwner::Vertex(vertex) => (feynkit_graph::symbols::vertex_index(), vertex),
            NumeratorOwner::Edge(edge) => (feynkit_graph::symbols::edge_index(), edge),
        };
        let owner = i64::try_from(owner)
            .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic numerator owner"))?;
        let local = i64::try_from(local)
            .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic dummy index"))?;
        if self.frame == NumeratorFrame::Interaction {
            Ok(self.dummy_symbol().call((0, local)))
        } else {
            Ok(symbol.call((owner, local)))
        }
    }

    fn momentum(
        &self,
        leg: NumeratorHalfEdge,
        index: Option<Atom>,
    ) -> Result<Atom, GenerationError> {
        let edge = i64::try_from(leg.edge)
            .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic edge identifier"))?;
        let mut momentum = FunctionBuilder::new(self.momentum_symbol()).add_arg(edge);
        if let Some(index) = index {
            momentum = momentum.add_arg(index);
        }
        Ok(momentum.finish())
    }

    /// Assign private rule slots to the finalized graph. Keep indices before
    /// carriers: replacing a carrier skips its arguments in a top-down walk.
    /// Half-edge slots stay distinct even when two momenta bind to a self-loop.
    fn bind_interaction(&self, expression: &Atom) -> Result<Atom, GenerationError> {
        let mut atom = expression.clone();
        let mut error = None;
        for carriers in [false, true] {
            atom = atom.replace_map(|term, _, out| {
                if error.is_some() {
                    return;
                }
                let AtomView::Fun(function) = term else {
                    return;
                };
                let head = function.get_symbol();
                if head.get_namespace() != "feynkit_generator_template" {
                    return;
                }
                let name = head.get_stripped_name();
                if matches!(name, "momentum" | "orientation") != carriers {
                    return;
                }
                let mut arguments = function.iter();
                let slot = usize::try_from(arguments.next().expect("compiled local slot"))
                    .expect("compiled local slot is nonnegative");
                let value = match name {
                    "hedge" => {
                        let shift = i64::try_from(arguments.next().expect("compiled index shift"))
                            .expect("compiled index shift is integral");
                        self.index(self.legs[slot], shift)
                    }
                    "dummy" => {
                        let local = usize::try_from(arguments.next().expect("compiled dummy"))
                            .expect("compiled dummy is nonnegative");
                        self.dummy(local)
                    }
                    "momentum" => self.momentum(
                        self.legs[slot],
                        arguments.next().map(|index| index.to_owned()),
                    ),
                    "orientation" => Ok(if self.legs[slot].flow == Flow::Source {
                        -Atom::one()
                    } else {
                        Atom::one()
                    }),
                    _ => unreachable!("private interaction template head"),
                };
                match value {
                    Ok(value) => **out = value,
                    Err(failure) => error = Some(failure),
                }
            });
        }
        error.map_or(Ok(atom), Err)
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

/// A rule's representation data is independent of graph positions and momentum flow.
#[derive(PartialEq, Eq, Hash)]
struct InteractionRuleKey {
    rule: VertexRuleId,
    colors: Vec<i64>,
}

/// Most UFO factors can be lowered before assigning graph indices. Slash dummy
/// allocation and dummy identity contraction must instead see the concrete
/// canonical expression order. Already graph-owned input has the same requirement.
enum InteractionFactor {
    Local(Atom),
    Deferred(Atom, NumeratorSector),
}

impl InteractionFactor {
    fn compile(
        expression: &Atom,
        sector: NumeratorSector,
        local: &NumeratorInstantiation<'_>,
    ) -> Self {
        let mut concrete = false;
        let mut identity = false;
        let mut dummy = false;
        expression.visitor(&mut |term| {
            dummy |= i64::try_from(term).is_ok_and(|index| index < 0);
            if let AtomView::Fun(function) = term {
                let head = function.get_symbol();
                let name = head.get_stripped_name();
                concrete |= [
                    feynkit_graph::symbols::hedge_index(),
                    feynkit_graph::symbols::vertex_index(),
                    feynkit_graph::symbols::edge_index(),
                    momentum_symbol(),
                ]
                .contains(&head)
                    || head.get_namespace() == "feynkit_generator_template";
                if NumeratorInstantiation::is_model_symbol(head, name) {
                    concrete |= name == "PSlash"
                        || (sector == NumeratorSector::Color && matches!(name, "P" | "Slash"));
                    identity |= name == "Identity";
                    dummy |= name == "dummy";
                }
            }
            true
        });
        concrete |= sector == NumeratorSector::Color && identity && dummy;
        if !concrete {
            // A malformed local factor can disappear after concrete momenta
            // coincide. Keep its validation at that concrete boundary too.
            if let Ok(atom) = local.instantiate(expression, sector) {
                return Self::Local(atom);
            }
        }
        Self::Deferred(expression.clone(), sector)
    }

    fn instantiate(&self, frame: &NumeratorInstantiation<'_>) -> Result<Atom, GenerationError> {
        match self {
            Self::Local(atom) => frame.bind_interaction(atom),
            Self::Deferred(atom, sector) => frame.instantiate(atom, *sector),
        }
    }
}

struct InteractionTerm {
    coupling: Atom,
    color: InteractionFactor,
    lorentz: InteractionFactor,
}

enum CompiledInteraction {
    Local(Atom),
    Mixed(Vec<InteractionTerm>),
}

impl CompiledInteraction {
    fn new(
        model: &Model,
        rule_id: VertexRuleId,
        frame: &NumeratorInstantiation<'_>,
    ) -> Result<Self, GenerationError> {
        let rule = model.vertex_rule_by_id(rule_id)?;
        let legs = frame
            .legs
            .iter()
            .enumerate()
            .map(|(slot, leg)| NumeratorHalfEdge {
                edge: slot,
                hedge: slot,
                ..*leg
            })
            .collect::<Vec<_>>();
        let local = NumeratorInstantiation {
            owner: frame.owner,
            legs: &legs,
            frame: NumeratorFrame::Interaction,
        };
        let mut terms = Vec::new();
        for (color_index, color) in rule.color_structures.iter().enumerate() {
            for (lorentz_index, lorentz_id) in rule.lorentz_structures.iter().enumerate() {
                let Some(coupling) = rule
                    .couplings
                    .get(color_index)
                    .and_then(|row| row.get(lorentz_index))
                    .and_then(Option::as_ref)
                else {
                    continue;
                };
                let coupling = model.coupling_by_id(*coupling)?;
                let coupling = Atom::var(symbol!(&format!("UFO::{}", coupling.name)));
                let lorentz = &model.lorentz_structure_by_id(*lorentz_id)?.structure;
                terms.push(InteractionTerm {
                    coupling,
                    color: InteractionFactor::compile(color, NumeratorSector::Color, &local),
                    lorentz: InteractionFactor::compile(lorentz, NumeratorSector::Spin, &local),
                });
            }
        }
        if terms.is_empty() {
            return Ok(Self::Local(Atom::one()));
        }
        if terms.iter().all(|term| {
            matches!(term.color, InteractionFactor::Local(_))
                && matches!(term.lorentz, InteractionFactor::Local(_))
        }) {
            Ok(Self::Local(Atom::add_many(terms.into_iter().map(|term| {
                let InteractionFactor::Local(color) = term.color else {
                    unreachable!()
                };
                let InteractionFactor::Local(lorentz) = term.lorentz else {
                    unreachable!()
                };
                term.coupling * color * lorentz
            }))))
        } else {
            Ok(Self::Mixed(terms))
        }
    }

    fn instantiate(&self, frame: &NumeratorInstantiation<'_>) -> Result<Atom, GenerationError> {
        match self {
            Self::Local(atom) => frame.bind_interaction(atom),
            Self::Mixed(terms) => {
                let mut instantiated = Vec::with_capacity(terms.len());
                for term in terms {
                    let color = term.color.instantiate(frame)?;
                    let lorentz = term.lorentz.instantiate(frame)?;
                    instantiated.push(&term.coupling * color * lorentz);
                }
                Ok(Atom::add_many(instantiated))
            }
        }
    }
}

/// Resolved process inputs and their one shared search builder. The search census and
/// projector oracle use the same pruning and budget construction as ordinary generation.
struct PreparedSearch {
    options: GenerationOptions,
    resolved: ResolvedProcess,
    signatures: VertexSignatures,
    externals: ProcessExternals,
    generator: GraphGenerator<ParticleId, ExternalLabel>,
}

#[derive(Clone)]
pub(crate) struct Generator {
    model: Arc<Model>,
    numerators: Arc<NumeratorCache>,
}

/// Localized factors depend on the immutable model and their index assignment,
/// not on the rest of the diagram. Keep them for this generator's lifetime.
#[derive(Default)]
struct NumeratorCache {
    // An exact frame hit also avoids binding the compiled rule again.
    interactions: RwLock<HashMap<InteractionNumeratorKey, Atom>>,
    interaction_rules: RwLock<HashMap<InteractionRuleKey, Arc<CompiledInteraction>>>,
    propagators: RwLock<HashMap<PropagatorNumeratorKey, Atom>>,
    /// Automatic projectors by the finalized position of every external state. Finalization
    /// numbers the external legs first, so all diagrams of an amplitude share one entry. A cross
    /// section adds the sewn final-state copies after the internal edges, so its entries vary
    /// with the internal edge count.
    projectors: RwLock<HashMap<Vec<ExternalStatePosition>, Atom>>,
}

/// The value cached under `key`, computed and stored on first use.
fn memoised<K: Eq + Hash, V: Clone>(
    cache: &RwLock<HashMap<K, V>>,
    key: K,
    compute: impl FnOnce(&K) -> Result<V, GenerationError>,
) -> Result<V, GenerationError> {
    if let Some(value) = cache.read().unwrap().get(&key) {
        return Ok(value.clone());
    }
    let value = compute(&key)?;
    cache.write().unwrap().insert(key, value.clone());
    Ok(value)
}

/// An external state at its position in a diagram: the edge that carries it and the half-edge
/// where it attaches.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct ExternalStatePosition {
    external: ExternalNode,
    edge: usize,
    hedge: usize,
}

impl ExternalStatePosition {
    /// The position of `external`, carried by the leg with `connection`, in a finalized diagram.
    /// An incoming state attaches at the target of its edge and an outgoing one at the source.
    /// `FeynmanDiagram::validate` requires this orientation of the unpaired legs. Finalization
    /// sews a cross-section initial state with the amplitude-side sink as the target and the
    /// conjugate-side source as the source, so the sewn edge carries both copies.
    fn finalized(external: ExternalNode, connection: usize, diagram: &FeynmanDiagram) -> Self {
        let (edge, _, _) = diagram
            .edges()
            .find(|(_, _, data)| {
                data.external
                    .as_ref()
                    .is_some_and(|leg| leg.connection == connection)
            })
            .expect("finalized external leg");
        let endpoint = match external.state {
            ExternalState::Incoming => DiagramEndpoint::Target,
            ExternalState::Outgoing => DiagramEndpoint::Source,
        };
        let hedge = diagram
            .half_edge_id(DiagramHalfEdge { edge, endpoint })
            .expect("finalized external endpoint");
        Self {
            external,
            edge: edge.0,
            hedge: hedge.0,
        }
    }
}

#[derive(PartialEq, Eq, Hash)]
struct InteractionNumeratorKey {
    rule: VertexRuleId,
    owner: NumeratorOwner,
    legs: Vec<NumeratorHalfEdge>,
}

#[derive(PartialEq, Eq, Hash)]
struct PropagatorNumeratorKey {
    particle: ParticleId,
    edge: usize,
    hedges: [usize; 2],
}

/// Feynman-rule numerators for a finalized topology, using its attached model.
pub trait FeynmanDiagramRulesExt {
    /// Instantiate vertex and propagator numerators in an imported diagram's
    /// existing half-edge frame, using its interaction-slot assignments.
    ///
    /// The model's Feynman rules replace local numerator fragments, just as during
    /// diagram generation, at every vertex with an interaction and on every internal
    /// edge. Dummy edges are no propagators or interaction legs and keep their
    /// fragments, as does a vertex without an interaction. The routing, overall
    /// factor and numerator prefactor are retained.
    ///
    /// A [`FeynmanDiagram`] keeps its projector and replaces all of these
    /// fragments. A [`DotDiagram`] keeps the fragments and the projector its input
    /// wrote and fills in the others: an omitted projector becomes the spinor or
    /// polarization vector of every external state, as in generation.
    fn apply_feynman_rules(&self) -> Result<FeynmanDiagram, GenerationError>;
}

impl FeynmanDiagramRulesExt for FeynmanDiagram {
    fn apply_feynman_rules(&self) -> Result<FeynmanDiagram, GenerationError> {
        // `map_data` returns unvalidated diagrams: check their slots before they
        // index the vertex rules. The check is cached for validated diagrams.
        self.validate()?;
        Generator::new(self.model_arc()).construct_numerators(self.clone(), |_| true)
    }
}

impl FeynmanDiagramRulesExt for DotDiagram {
    fn apply_feynman_rules(&self) -> Result<FeynmanDiagram, GenerationError> {
        self.diagram.validate()?;
        let generator = Generator::new(self.diagram.model_arc());
        let diagram =
            generator.construct_numerators(self.diagram.clone(), |owner| match owner {
                NumeratorOwner::Vertex(vertex) => {
                    !self.explicit.vertices.contains(&VertexId(vertex))
                }
                NumeratorOwner::Edge(edge) => !self.explicit.edges.contains(&EdgeId(edge)),
            })?;
        if self.explicit.projector {
            return Ok(diagram);
        }
        let projector = generator.finalized_external_state_projector(&diagram)?;
        Ok(diagram.with_projector(projector))
    }
}

impl Generator {
    pub fn new(model: impl Into<Arc<Model>>) -> Self {
        Self {
            model: model.into(),
            numerators: Arc::default(),
        }
    }

    fn search(
        &self,
        process: &Process,
        options: &GenerationOptions,
        generation_type: GenerationType,
        observe: impl Fn(&GenerationView<'_>, bool) + Send + Sync + 'static,
    ) -> Result<PreparedSearch, GenerationError> {
        process.validate()?;
        if generation_type == GenerationType::Amplitude
            && process.outgoing_alternatives().len() != 1
        {
            return Err(ProcessError::MultipleAmplitudeFinalStates.into());
        }
        process.validate_in(&self.model)?;
        let options = process
            .restrict_options(&self.model, options)?
            .resolve_selectors(&self.model)?;
        self.validate_options(generation_type, &options)?;
        process.validate_covariant_cut_filters(&self.model, &options, generation_type)?;
        let expanded_process = process.clone().with_final_state_alternatives(
            process.covariant_cut_states(&self.model, generation_type)?,
        )?;
        let process = &expanded_process;
        let resolved = ResolvedProcess::new(&self.model, process, generation_type, &options)?;
        let signatures = VertexSignatures::new(&self.model, &options)?;
        let externals = resolved.externals(&self.model)?;
        let pruner = Pruner::new(&self.model, &options);

        let generator = signatures.generator(&self.model, &externals, &options, move |view| {
            let kept = pruner.keeps(view);
            observe(view, kept);
            kept
        })?;
        Ok(PreparedSearch {
            options,
            resolved,
            signatures,
            externals,
            generator,
        })
    }

    pub fn generate(
        &self,
        process: &Process,
        options: &GenerationOptions,
        generation_type: GenerationType,
    ) -> Result<GenerationResult, GenerationError> {
        let PreparedSearch {
            options,
            resolved,
            signatures,
            externals,
            generator,
        } = self.search(process, options, generation_type, |_, _| {})?;
        let options = &options;

        options.report_progress("topologies", 0, None);
        let topologies = Arc::new(AtomicUsize::new(0));
        let generator = generator
            .progress({
                let options = options.clone();
                let topologies = Arc::clone(&topologies);
                move |progress| {
                    topologies.store(progress.topologies, Ordering::Relaxed);
                    options.report_progress("topologies", progress.topologies, None);
                    !options.cancellation_requested()
                }
            })
            .abort_check({
                let options = options.clone();
                move || options.cancellation_requested()
            });
        #[cfg(not(target_arch = "wasm32"))]
        let generator = generator.parallel(true);
        let enumerate = || {
            generator
                .iter()
                .map(|diagrams| diagrams.collect::<Vec<_>>())
        };
        #[cfg(not(target_arch = "wasm32"))]
        let pool = {
            let mut pool = ThreadPoolBuilder::new();
            // Respect an explicit thread count; otherwise Rayon uses its default
            // logical-CPU count. Browser kernels process the same graphs serially.
            if let Some(threads) = options.threads {
                pool = pool.num_threads(threads);
            }
            pool.build()?
        };
        let mut generated = if externals.vetoed(&self.model, options) {
            Vec::new()
        } else {
            #[cfg(target_arch = "wasm32")]
            let generated = enumerate()?;
            #[cfg(not(target_arch = "wasm32"))]
            let generated = pool.install(enumerate)?;
            generated
        };
        let mut completed = !options.cancellation_requested();
        // Numerator grouping keeps the first representative, so preserve its canonical
        // encounter order. Ungrouped output retains Linnet's deterministic search order
        // without demanding canonical forms merely for naming.
        if !matches!(options.numerator_grouping, NumeratorGrouping::None) {
            generated.sort_by(|a, b| a.canonical_form().cmp(b.canonical_form()));
        }
        let topology_count = topologies.load(Ordering::Relaxed);
        let interaction_assignment_count = generated.len();

        let prepared_count = Mutex::new(0);
        options.report_progress("interaction_filters", 0, Some(interaction_assignment_count));
        let prepare = || {
            #[cfg(target_arch = "wasm32")]
            let diagrams = generated.iter();
            #[cfg(not(target_arch = "wasm32"))]
            let diagrams = generated.par_iter();
            diagrams
                .map(|diagram| {
                    if options.cancellation_requested() {
                        return Ok(None);
                    }
                    // Diagrams related by permuting sewn pairs are kept once, weighted by
                    // their number.
                    let input = match externals.sewn_pair_orbit(diagram) {
                        None => None,
                        Some(multiplicity) => {
                            let graph = externals.graph(diagram, &signatures);
                            let symmetry = diagram
                                .automorphism_count
                                .and_then(|count| u64::try_from(count).ok())
                                .ok_or(GenerationError::SymmetryFactor)?;
                            self.diagram_input(&resolved, graph, symmetry, multiplicity, options)?
                        }
                    };
                    let mut count = prepared_count.lock().unwrap();
                    *count += 1;
                    options.report_progress(
                        "interaction_filters",
                        *count,
                        Some(interaction_assignment_count),
                    );
                    Ok(input)
                })
                .collect::<Result<Vec<_>, GenerationError>>()
        };
        #[cfg(target_arch = "wasm32")]
        let diagram_inputs = prepare()?;
        #[cfg(not(target_arch = "wasm32"))]
        let diagram_inputs = pool.install(prepare)?;
        let mut diagram_inputs = diagram_inputs.into_iter().flatten().collect::<Vec<_>>();
        if resolved.identifies_mirror_images() {
            diagram_inputs = DiagramInput::merge_identical(diagram_inputs);
        }
        if options.cancellation_requested() {
            completed = false;
        }
        let total = diagram_inputs.len();
        options.report_progress("numerators", 0, Some(total));
        let constructed = Mutex::new(0);
        let construct = || {
            #[cfg(target_arch = "wasm32")]
            let inputs = diagram_inputs.into_iter();
            #[cfg(not(target_arch = "wasm32"))]
            let inputs = diagram_inputs.into_par_iter();
            inputs
                .map(|input| {
                    if options.cancellation_requested() {
                        return Ok::<_, GenerationError>(None);
                    }
                    let pair = self.diagram_pair(None, generation_type, input, options)?;
                    let mut count = constructed.lock().unwrap();
                    *count += 1;
                    options.report_progress("numerators", *count, Some(total));
                    Ok(Some(pair))
                })
                .collect::<Result<Vec<_>, _>>()
        };
        #[cfg(target_arch = "wasm32")]
        let diagram_pairs = construct()?;
        #[cfg(not(target_arch = "wasm32"))]
        let diagram_pairs = pool.install(construct)?;
        let mut diagram_pairs = diagram_pairs.into_iter().flatten().collect::<Vec<_>>();
        if options.cancellation_requested() {
            completed = false;
        }
        options.report_progress("selection", 0, None);
        if !options.forced_cuts.is_empty() {
            let mut matched = vec![false; options.forced_cuts.len()];
            let mut retained = Vec::new();
            for pair in diagram_pairs {
                let cuts = pair
                    .comparison()
                    .cuts()
                    .iter()
                    .filter(|cut| {
                        let edges = cut
                            .cut
                            .iter()
                            .map(|half_edge| half_edge.edge)
                            .collect::<BTreeSet<_>>();
                        let keep = options
                            .forced_cuts
                            .iter()
                            .enumerate()
                            .filter(|(_, requested)| **requested == edges)
                            .map(|(index, _)| index)
                            .collect::<Vec<_>>();
                        for index in keep {
                            matched[index] = true;
                        }
                        options.forced_cuts.contains(&edges)
                    })
                    .cloned()
                    .collect::<Vec<_>>();
                if !cuts.is_empty() {
                    let pair = pair.with_cuts(cuts)?;
                    retained.push(pair);
                }
            }
            diagram_pairs = retained;
            if let Some((_, missing)) = options
                .forced_cuts
                .iter()
                .enumerate()
                .find(|(index, _)| !matched[*index])
            {
                return Err(GenerationError::UnknownForcedCut {
                    edges: missing.iter().map(|edge| edge.0).collect(),
                });
            }
        }
        let naming = match &options.output {
            GenerationOutput::Minimal => None,
            GenerationOutput::Named { prefix } => Some((
                prefix,
                diagram_pairs.len().saturating_sub(1).to_string().len(),
            )),
        };
        let mut selected = Vec::with_capacity(diagram_pairs.len());
        for (index, mut pair) in diagram_pairs.into_iter().enumerate() {
            let name = naming.map(|(prefix, width)| format!("{prefix}{index:0width$}"));
            let selected_by_id = match &options.selected_diagram_ids {
                Some(selected) if !selected.is_empty() => {
                    selected.contains(&pair.comparison().id()?)
                        || selected.contains(&pair.representative().id()?)
                }
                _ => false,
            };
            let selected_by_name = name.as_ref().is_some_and(|name| {
                options
                    .selected_diagram_names
                    .as_ref()
                    .is_some_and(|selected| selected.contains(name))
            });
            // A present empty selection still selects nothing, without demanding an ID.
            let has_selection =
                options.selected_diagram_ids.is_some() || options.selected_diagram_names.is_some();
            if (has_selection && !selected_by_id && !selected_by_name)
                || (!options.vetoed_diagram_ids.is_empty()
                    && (options
                        .vetoed_diagram_ids
                        .contains(&pair.comparison().id()?)
                        || options
                            .vetoed_diagram_ids
                            .contains(&pair.representative().id()?)))
                || name
                    .as_ref()
                    .is_some_and(|name| options.vetoed_diagram_names.contains(name))
            {
                continue;
            }
            if let Some(name) = &name {
                pair = pair.with_name(name);
            }
            let by_id = if options.loop_momentum_bases.is_empty() {
                None
            } else if let Some(edges) = options
                .loop_momentum_bases
                .get(&pair.representative().id()?)
            {
                Some(edges)
            } else {
                options.loop_momentum_bases.get(&pair.comparison().id()?)
            };
            let requested_basis = by_id.or_else(|| {
                name.as_ref()
                    .and_then(|name| options.named_loop_momentum_bases.get(name))
            });
            if let Some(edges) = requested_basis {
                pair = pair.with_loop_momentum_edges(edges)?;
            }
            selected.push(pair);
        }
        let grouping = if options.cancellation_requested() {
            completed = false;
            &NumeratorGrouping::None
        } else {
            &options.numerator_grouping
        };
        let group = || {
            grouping::group_diagrams(
                selected,
                &self.model,
                grouping,
                options.symmetrizes_left_right(),
                options,
            )
        };
        #[cfg(target_arch = "wasm32")]
        let grouped = group()?;
        #[cfg(not(target_arch = "wasm32"))]
        let grouped = pool.install(group)?;
        if options.cancellation_requested() {
            completed = false;
        }

        let zero_numerator_count = grouped.zero_numerator_count;
        let groups = grouped.groups;
        let diagrams = grouped.diagrams;
        let total = groups.len();
        options.report_progress_with_zero_count("grouping", 0, Some(total), zero_numerator_count);
        for output_index in 0..total {
            options.report_progress_with_zero_count(
                "grouping",
                output_index + 1,
                Some(total),
                zero_numerator_count,
            );
        }

        let result = GenerationResult {
            report: GenerationReport {
                topology_count,
                interaction_assignment_count,
                retained_count: diagrams.len(),
                zero_numerator_count,
                completed,
            },
            diagrams,
            groups,
        };
        result.validate_groups()?;
        options.report_progress_with_zero_count(
            if completed { "complete" } else { "cancelled" },
            result.diagrams.len(),
            Some(result.diagrams.len()),
            zero_numerator_count,
        );
        Ok(result)
    }

    /// Canonicalize a generated diagram, apply the cut and interaction filters that need its
    /// canonical frame, and normalize its fermion flows.
    fn diagram_input(
        &self,
        resolved: &ResolvedProcess,
        graph: GeneratedGraph,
        symmetry: u64,
        multiplicity: u64,
        options: &GenerationOptions,
    ) -> Result<Option<DiagramInput>, GenerationError> {
        let mut topology = ColoredTopology {
            graph: resolved.canonicalize_numerator_graph(
                &graph,
                &self.model,
                &options.numerator_grouping,
            )?,
            symmetry,
            cut_partitions: Vec::new(),
            topology_threshold_candidates: Vec::new(),
        };
        if resolved.generation_type == GenerationType::CrossSection {
            let inventory = resolved.cut_partitions(&topology.graph, &self.model, options)?;
            topology.cut_partitions = inventory.physical;
            topology.topology_threshold_candidates = inventory.topology_threshold_candidates;
            if topology.cut_partitions.is_empty()
                || !self.retain_cut_amplitude_partitions(&mut topology, options)?
                || !resolved.passes_sewn_filter(&topology.graph, options)
            {
                return Ok(None);
            }
        }
        topology.graph = ResolvedProcess::repair_cp_interaction_rules(topology.graph, &self.model)?;
        if !self.passes_vertex_selection(&topology.graph, options) {
            return Ok(None);
        }
        let grassmann_loop_count = self.closed_fermion_loop_count(&topology.graph, true)?;
        // Canonicalization fixes the stable vertex/edge-ID frame used by
        // cuts and routing.  The physical representative must still have
        // each fermion chain oriented consistently before its numerator
        // is instantiated.  Reversing an edge keeps its stable ID, so the
        // only cut metadata that needs transport is the source/target
        // endpoint attached to each half-edge.
        let normalized = resolved.normalize_fermion_flows(&topology.graph, &self.model)?;
        let representative = ColoredTopology {
            graph: normalized.graph,
            symmetry,
            cut_partitions: transport_cut_metadata(
                &topology.cut_partitions,
                &normalized.reversed_edges,
            ),
            topology_threshold_candidates: transport_threshold_candidate_metadata(
                &topology.topology_threshold_candidates,
                &normalized.reversed_edges,
            ),
        };
        // Without numerator comparisons, the normalized physical graph is
        // already the representative.
        let (comparison, reversed_edges) =
            if matches!(options.numerator_grouping, NumeratorGrouping::None) {
                (representative.clone(), BTreeSet::new())
            } else {
                (topology, normalized.reversed_edges)
            };
        Ok(Some(DiagramInput {
            comparison,
            representative,
            signs: FermionSigns {
                internal_loop_negative: grassmann_loop_count % 2 == 1,
                ..normalized.signs
            },
            reversed_edges,
            multiplicity,
        }))
    }

    /// Finalize the comparison and representative diagrams of a prepared input.
    fn diagram_pair(
        &self,
        name: Option<String>,
        generation_type: GenerationType,
        input: DiagramInput,
        options: &GenerationOptions,
    ) -> Result<DiagramPair, GenerationError> {
        let DiagramInput {
            comparison,
            representative,
            signs,
            reversed_edges,
            multiplicity,
        } = input;
        let identical = comparison == representative;
        let internal_reversals = comparison
            .graph
            .edges()
            .enumerate()
            .filter(|(_, edge)| !edge.is_external())
            .map(|(edge, _)| reversed_edges.contains(&EdgeId(edge)))
            .collect::<Vec<_>>();
        let comparison = self.to_diagram(
            name.clone(),
            generation_type,
            comparison,
            signs,
            multiplicity,
            options,
        )?;
        if identical {
            return Ok(DiagramPair::Shared(comparison));
        }
        let representative = self.to_diagram(
            name,
            generation_type,
            representative,
            signs,
            multiplicity,
            options,
        )?;
        // Finalization orders external carriers before the internal edges, which keep the
        // canonical order. Transport the comparison-to-representative reversals into that
        // finalized edge frame; process flow fixes all external carriers.
        let reversed_edges = comparison
            .edges()
            .filter(|(_, _, edge)| edge.external.is_none())
            .zip(internal_reversals)
            .filter_map(|((edge, _, _), reversed)| reversed.then_some(edge))
            .collect();
        Ok(DiagramPair::Distinct {
            comparison,
            representative: Box::new(representative),
            reversed_edges,
        })
    }

    fn validate_options(
        &self,
        generation_type: GenerationType,
        options: &GenerationOptions,
    ) -> Result<(), GenerationError> {
        if options.threads == Some(0) {
            return Err(GenerationError::InvalidThreadCount(0));
        }
        options.validate_for(generation_type)?;
        if matches!(
            &options.numerator_grouping,
            NumeratorGrouping::Identical(options)
                | NumeratorGrouping::UpToSign(options)
                | NumeratorGrouping::UpToScalar(options)
                if options.number_of_numerical_samples == 0
        ) {
            return Err(GenerationError::InvalidNumericalSampleCount);
        }
        self.validate_filters(FilterScope::Graph, generation_type, &options.graph_filters)?;
        self.validate_filters(
            FilterScope::CutAmplitude,
            generation_type,
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
                    | GenerationFilter::LoopOneParticleIrreducible
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
                        | GenerationFilter::LoopOneParticleIrreducible
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
                GenerationFilter::ParticleVeto(selectors) => {
                    for selector in selectors {
                        selector.resolve(&self.model)?;
                    }
                }
                GenerationFilter::VertexAllow(selectors)
                | GenerationFilter::VertexVeto(selectors) => {
                    for selector in selectors {
                        selector.resolve(&self.model)?;
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
                | GenerationFilter::LoopOneParticleIrreducible
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
            | GenerationFilter::LoopOneParticleIrreducible
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

    /// The vertex selections, rechecked after charge-conjugation repair: the generator only
    /// offers selected rules, but a repaired rule may be excluded.
    fn passes_vertex_selection(&self, graph: &GeneratedGraph, options: &GenerationOptions) -> bool {
        options.graph_filters.iter().all(|filter| match filter {
            GenerationFilter::VertexAllow(selectors) => graph.vertices().all(|rule| {
                selectors
                    .iter()
                    .any(|selector| resolved_vertex(selector) == rule)
            }),
            GenerationFilter::VertexVeto(selectors) => graph.vertices().all(|rule| {
                !selectors
                    .iter()
                    .any(|selector| resolved_vertex(selector) == rule)
            }),
            _ => true,
        })
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

    // Count closed, non-branching fermion components, optionally including ghosts
    // for the Grassmann sign. The fermion-loop filter excludes ghosts. Branching
    // components are rejected explicitly.
    fn closed_fermion_loop_count(
        &self,
        graph: &GeneratedGraph,
        include_ghosts: bool,
    ) -> Result<usize, GenerationError> {
        let nodes = graph.plain_nodes();
        let mut ends = Vec::new();
        for (id, edge) in graph.edges().enumerate() {
            let particle = self.model.particle_by_id(edge.particle)?;
            if particle.is_fermion() && !particle.is_ghost()
                || include_ghosts && particle.is_ghost()
            {
                ends.push(nodes.ends(id, &edge));
            }
        }
        closed_loops(nodes.len(), &ends).map_err(|(vertex, degree)| {
            GenerationError::UnsupportedFermionBranching { vertex, degree }
        })
    }

    fn to_diagram(
        &self,
        name: Option<String>,
        generation_type: GenerationType,
        topology: ColoredTopology,
        fermion_signs: FermionSigns,
        multiplicity: u64,
        options: &GenerationOptions,
    ) -> Result<FeynmanDiagram, GenerationError> {
        // The automorphism count of the rule-coloured graph is a denominator
        // in the overall multiplier.
        let mut factors = vec![
            symbol!("feynkit_generator_factor::AutG")
                .call(topology.symmetry)
                .pow(-1),
        ];
        if multiplicity > 1 {
            factors.push(
                symbol!("feynkit_generator_factor::NumeratorIndependentSymmetryGrouping")
                    .call(multiplicity),
            );
        }
        if fermion_signs.internal_loop_negative {
            factors.push(symbol!("feynkit_generator_factor::InternalFermionLoopSign").call(-1));
        }
        if fermion_signs.include_external_ordering {
            let sign = if fermion_signs.external_ordering_negative {
                -1
            } else {
                1
            };
            factors
                .push(symbol!("feynkit_generator_factor::ExternalFermionOrderingSign").call(sign));
        }
        if fermion_signs.include_antifermion_spin_sum {
            let sign = if fermion_signs.antifermion_spin_sum_negative {
                -1
            } else {
                1
            };
            factors.push(symbol!("feynkit_generator_factor::AntiFermionSpinSumSign").call(sign));
        }
        let overall_factor = Atom::mul_many(factors);
        let graph = &topology.graph;
        let incoming_count = graph
            .externals()
            .filter(|(_, external, _)| external.state == ExternalState::Incoming)
            .count();
        let mut builder = FeynmanDiagramBuilder::new(Arc::clone(&self.model), name)
            .canonical_loop_momentum_basis()
            .symmetry_factor(topology.symmetry)
            .overall_factor(overall_factor)
            .numerator_prefactor(options.numerator_prefactor.clone())
            .cuts(topology.cut_partitions.clone())
            .topology_threshold_candidates(topology.topology_threshold_candidates.clone());
        // Every external state with the connection of the leg that carries it.
        let mut external_states = Vec::new();
        let mut endpoint_slots = vec![[VertexSlot(0); 2]; graph.n_edges()];
        for (vertex, rule_id) in graph.vertices().enumerate() {
            let rule = self.model.vertex_rule_by_id(rule_id)?;
            let legs = self.interaction_half_edges(graph, vertex, rule)?;
            for (slot, leg) in legs.iter().enumerate() {
                endpoint_slots[leg.edge][usize::from(leg.flow.is_sink())] = VertexSlot(slot);
            }
            builder.add_vertex(DiagramVertex::interaction(format!("v{vertex}"), rule_id));
        }
        for (edge_index, edge) in graph.edges().enumerate() {
            let diagram_edge = DiagramEdge::new(edge.particle, edge.directed);
            let slot = |flow: Flow| endpoint_slots[edge_index][usize::from(flow.is_sink())];
            match (edge.source, edge.sink) {
                (Endpoint::Vertex(source), Endpoint::Vertex(sink)) => {
                    builder.add_edge_with_slots(
                        feynkit_graph::VertexId(source),
                        feynkit_graph::VertexId(sink),
                        diagram_edge,
                        slot(Flow::Source),
                        slot(Flow::Sink),
                    )?;
                }
                (Endpoint::External(external), Endpoint::Vertex(vertex))
                | (Endpoint::Vertex(vertex), Endpoint::External(external)) => {
                    let flow = if edge.sink == Endpoint::Vertex(vertex) {
                        Flow::Sink
                    } else {
                        Flow::Source
                    };
                    let connection = match (generation_type, external.state) {
                        (GenerationType::Amplitude, _) | (_, ExternalState::Incoming) => {
                            external.index
                        }
                        (GenerationType::CrossSection, ExternalState::Outgoing) => external
                            .index
                            .checked_sub(incoming_count)
                            .ok_or(GenerationError::MissingExternalParticle {
                                index: external.index,
                            })?,
                    };
                    external_states.push((external, connection));
                    builder.add_generation_external(
                        feynkit_graph::VertexId(vertex),
                        flow,
                        diagram_edge,
                        slot(flow),
                        ExternalLeg {
                            name: format!("ext{}", external.index),
                            index: external.index,
                            state: external.state,
                            connection,
                        },
                    )?;
                }
                (Endpoint::External(_), Endpoint::External(_)) => {
                    unreachable!("generated edges touch an interaction vertex")
                }
            }
        }
        let diagram = match &options.projector {
            // A requested projector refers to the provisional frame and is finalized with the
            // diagram.
            Some(projector) => builder.projector(projector.clone()).build()?,
            // The automatic one is built directly in the final frame.
            None => {
                let diagram = builder.build()?;
                let positions = external_states
                    .into_iter()
                    .map(|(external, connection)| {
                        ExternalStatePosition::finalized(external, connection, &diagram)
                    })
                    .collect();
                let projector = memoised(&self.numerators.projectors, positions, |positions| {
                    self.external_state_projector(positions)
                })?;
                diagram.with_projector(projector)
            }
        };
        self.construct_numerators(diagram, |_| true)
    }

    /// Instantiate the Feynman rules of the fragments `replaced` selects.
    fn construct_numerators(
        &self,
        diagram: FeynmanDiagram,
        replaced: impl Fn(NumeratorOwner) -> bool,
    ) -> Result<FeynmanDiagram, GenerationError> {
        // Instantiate directly in the finalized frame. This avoids building
        // every factor in provisional indices and then translating it again.
        let mut vertices = Vec::new();
        let graph = diagram.underlying();
        for (vertex, incident, data) in graph.iter_nodes() {
            let Some(rule_id) = data.interaction else {
                continue;
            };
            if !replaced(NumeratorOwner::Vertex(vertex.0)) {
                continue;
            }
            let rule = self.model.vertex_rule_by_id(rule_id)?;
            let mut legs = vec![None; rule.particles.len()];
            // Finalization carries incidence and UFO slots into the native half-edge frame.
            // Visit both ends of a self-loop, skipping dummy slots after the interaction legs.
            for hedge in incident {
                let edge = graph[&hedge];
                let data = &graph[edge];
                if data.is_dummy {
                    continue;
                }
                let flow = graph.flow(hedge);
                let slot = match flow {
                    Flow::Source => data.source_slot(),
                    Flow::Sink => data.target_slot(),
                };
                let particle = self.model.particle_by_id(rule.particles[slot.0])?;
                legs[slot.0] = Some(NumeratorHalfEdge {
                    edge: edge.0,
                    hedge: hedge.0,
                    flow,
                    spin: particle.spin,
                    color: particle.color,
                });
            }
            let legs = legs
                .into_iter()
                .map(|leg| leg.expect("assigned interaction slot"))
                .collect::<Vec<_>>();
            let instantiation =
                NumeratorInstantiation::new(NumeratorOwner::Vertex(vertex.0), &legs);
            vertices.push((
                VertexId(vertex.0),
                self.interaction_numerator(rule_id, &instantiation)?,
            ));
        }
        let mut edges = Vec::new();
        for (edge, _, data) in diagram.edges().filter(|(edge, _, data)| {
            data.external.is_none() && !data.is_dummy && replaced(NumeratorOwner::Edge(edge.0))
        }) {
            let hedges = [DiagramEndpoint::Source, DiagramEndpoint::Target].map(|endpoint| {
                diagram
                    .half_edge_id(DiagramHalfEdge { edge, endpoint })
                    .expect("internal edge endpoint")
                    .0
            });
            edges.push((edge, self.edge_numerator(data.particle, edge.0, hedges)?));
        }
        Ok(diagram.with_numerator_fragments(vertices, edges)?)
    }

    /// The external-state projector: the spinor or polarization vector of every state, on the
    /// edge that carries it and indexed at the half-edge where it attaches. Scalars contribute no
    /// factor.
    fn external_state_projector(
        &self,
        positions: &[ExternalStatePosition],
    ) -> Result<Atom, GenerationError> {
        let mut projector = Atom::one();
        for &ExternalStatePosition {
            external,
            edge,
            hedge,
        } in positions
        {
            if let Some(factor) =
                self.external_state_factor(external.particle, external.state, edge, hedge)?
            {
                projector *= factor;
            }
        }
        Ok(projector)
    }

    /// The external-state projector of a finalized diagram. An incoming state is
    /// indexed at the target of its edge and an outgoing one at the source; a sewn
    /// cross-section initial state carries both copies.
    fn finalized_external_state_projector(
        &self,
        diagram: &FeynmanDiagram,
    ) -> Result<Atom, GenerationError> {
        let mut projector = Atom::one();
        for (edge, endpoints, data) in diagram.edges() {
            let Some(leg) = &data.external else {
                continue;
            };
            let states = if endpoints.source.is_some() && endpoints.target.is_some() {
                vec![ExternalState::Incoming, ExternalState::Outgoing]
            } else {
                vec![leg.state]
            };
            for state in states {
                let endpoint = match state {
                    ExternalState::Incoming => DiagramEndpoint::Target,
                    ExternalState::Outgoing => DiagramEndpoint::Source,
                };
                let hedge = diagram
                    .half_edge_id(DiagramHalfEdge { edge, endpoint })
                    .expect("finalized external endpoint");
                if let Some(factor) =
                    self.external_state_factor(data.particle, state, edge.0, hedge.0)?
                {
                    projector *= factor;
                }
            }
        }
        Ok(projector)
    }

    /// The spinor or polarization vector of one external state on `edge`, indexed at
    /// `hedge`. Scalars have none.
    fn external_state_factor(
        &self,
        particle: ParticleId,
        state: ExternalState,
        edge: usize,
        hedge: usize,
    ) -> Result<Option<Atom>, GenerationError> {
        let index = feynkit_graph::symbols::hedge_index().call((hedge, 1_i64));
        let particle = self.model.particle_by_id(particle)?;
        let (head, index) = match particle.spin {
            2 => {
                let head = match (state, particle.is_antiparticle()) {
                    (ExternalState::Incoming, false) => feynkit_graph::symbols::u(),
                    (ExternalState::Incoming, true) => feynkit_graph::symbols::vbar(),
                    (ExternalState::Outgoing, false) => feynkit_graph::symbols::ubar(),
                    (ExternalState::Outgoing, true) => feynkit_graph::symbols::v(),
                };
                (head, Bispinor {}.new_rep(4).to_symbolic([index]))
            }
            3 => {
                let head = match state {
                    ExternalState::Incoming => feynkit_graph::symbols::epsilon(),
                    ExternalState::Outgoing => feynkit_graph::symbols::epsilonbar(),
                };
                (head, Minkowski {}.new_rep(4).to_symbolic([index]))
            }
            _ => return Ok(None),
        };
        Ok(Some(
            FunctionBuilder::new(head)
                .add_arg(edge)
                .add_arg(index)
                .finish(),
        ))
    }

    fn interaction_half_edges(
        &self,
        graph: &GeneratedGraph,
        vertex: usize,
        rule: &VertexRule,
    ) -> Result<Vec<NumeratorHalfEdge>, GenerationError> {
        // Reconstruct the UFO particle order discarded by topology
        // canonicalization. Direction distinguishes particle and antiparticle
        // slots, while a self-loop contributes both of its half-edges.
        let mut available = available_interaction_half_edges(graph, vertex, &self.model)?;

        let mut ordered = Vec::with_capacity(rule.particles.len());
        for (leg, particle_id) in rule.particles.iter().enumerate() {
            let particle = self.model.particle_by_id(*particle_id)?;
            let expected = incoming_stub(&self.model, particle);
            let Some(position) = available
                .iter()
                .position(|(_, signature)| *signature == expected)
            else {
                return Err(GenerationError::MissingInteractionLeg {
                    vertex,
                    interaction: rule.name.clone(),
                    leg: leg + 1,
                    particle: particle.name.clone(),
                });
            };
            let mut selected = available.remove(position).0;
            selected.spin = particle.spin;
            selected.color = particle.color;
            ordered.push(selected);
        }
        Ok(ordered)
    }

    fn propagator_numerator(&self, particle: &Particle) -> Result<(Atom, bool), GenerationError> {
        match self.model.particle_propagator(particle)? {
            Some((id, reversed)) => {
                Ok((self.model.propagator_by_id(id)?.numerator.clone(), reversed))
            }
            None => Ok((Atom::one(), false)),
        }
    }

    fn interaction_numerator(
        &self,
        rule_id: VertexRuleId,
        instantiation: &NumeratorInstantiation<'_>,
    ) -> Result<Atom, GenerationError> {
        let key = InteractionNumeratorKey {
            rule: rule_id,
            owner: instantiation.owner,
            legs: instantiation.legs.to_vec(),
        };
        memoised(&self.numerators.interactions, key, |_| {
            let rule_key = InteractionRuleKey {
                rule: rule_id,
                colors: instantiation.legs.iter().map(|leg| leg.color).collect(),
            };
            let rule = memoised(&self.numerators.interaction_rules, rule_key, |_| {
                CompiledInteraction::new(&self.model, rule_id, instantiation).map(Arc::new)
            })?;
            rule.instantiate(instantiation)
        })
    }

    fn edge_numerator(
        &self,
        particle_id: ParticleId,
        edge: usize,
        hedges: [usize; 2],
    ) -> Result<Atom, GenerationError> {
        let key = PropagatorNumeratorKey {
            particle: particle_id,
            edge,
            hedges,
        };
        memoised(&self.numerators.propagators, key, |_| {
            let particle = self.model.particle_by_id(particle_id)?;
            let (propagator, reversed) = self.propagator_numerator(particle)?;
            let flows = if reversed {
                [Flow::Sink, Flow::Source]
            } else {
                [Flow::Source, Flow::Sink]
            };
            let legs = flows.map(|flow| NumeratorHalfEdge {
                edge,
                flow,
                hedge: hedges[usize::from(flow == Flow::Sink)],
                spin: particle.spin,
                color: particle.color,
            });
            let instantiation = NumeratorInstantiation::new(NumeratorOwner::Edge(edge), &legs);
            let spin = instantiation.instantiate(&propagator, NumeratorSector::Spin)?;
            Ok(match instantiation.edge_color_identity(particle.color)? {
                None => spin,
                Some(color) => spin * color,
            })
        })
    }
}

/// The incident half-edges of `vertex` with the stub each one presents to the vertex: the
/// particle entering the vertex through it, as in the all-incoming UFO convention.
fn available_interaction_half_edges(
    graph: &GeneratedGraph,
    vertex: usize,
    model: &Model,
) -> Result<Vec<(NumeratorHalfEdge, ParticleStub)>, GenerationError> {
    graph
        .incident(vertex)
        .into_iter()
        .map(|(edge_id, flow)| {
            let edge = graph.edge(edge_id);
            let particle = model.particle_by_id(edge.particle)?;
            let entering = if flow == Flow::Sink || !edge.directed {
                particle
            } else {
                model.antiparticle(particle)?
            };
            Ok((
                NumeratorHalfEdge {
                    edge: edge_id,
                    flow,
                    hedge: 2 * edge_id + usize::from(flow == Flow::Sink),
                    spin: 0,
                    color: 1,
                },
                incoming_stub(model, entering),
            ))
        })
        .collect()
}

fn vertex_rule_signature(
    model: &Model,
    rule: &VertexRule,
) -> Result<Vec<ParticleStub>, GenerationError> {
    let mut signature = rule
        .particles
        .iter()
        .map(|particle| Ok(incoming_stub(model, model.particle_by_id(*particle)?)))
        .collect::<Result<Vec<_>, GenerationError>>()?;
    signature.sort();
    Ok(signature)
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct ColoredTopology {
    graph: GeneratedGraph,
    symmetry: u64,
    cut_partitions: Vec<DiagramCut>,
    topology_threshold_candidates: Vec<DiagramThresholdCandidate>,
}

/// A canonicalized diagram ready for finalization: the frame used for numerator comparisons,
/// the physical representative with normalized fermion flows, and the edges reversed between
/// them.
struct DiagramInput {
    comparison: ColoredTopology,
    representative: ColoredTopology,
    signs: FermionSigns,
    reversed_edges: BTreeSet<EdgeId>,
    /// The number of exactly labelled generated diagrams the representative stands for.
    multiplicity: u64,
}

impl DiagramInput {
    /// Merge the inputs with equal representatives and fermion signs into the first of them,
    /// adding up their multiplicities.
    ///
    /// Left/right symmetrization canonicalizes a sewn graph and its mirror image to the same
    /// representative. Kept apart, they would be two diagrams with one [`DiagramId`].
    fn merge_identical(inputs: Vec<Self>) -> Vec<Self> {
        let mut firsts = HashMap::with_capacity(inputs.len());
        let first = inputs
            .iter()
            .enumerate()
            .map(|(index, input)| {
                *firsts
                    .entry((&input.representative, input.signs))
                    .or_insert(index)
            })
            .collect::<Vec<_>>();
        let mut multiplicity = vec![0; inputs.len()];
        for (input, &first) in inputs.iter().zip(&first) {
            multiplicity[first] += input.multiplicity;
        }
        inputs
            .into_iter()
            .enumerate()
            .filter(|&(index, _)| first[index] == index)
            .map(|(index, input)| Self {
                multiplicity: multiplicity[index],
                ..input
            })
            .collect()
    }
}

/// Equal comparison and physical frames have one owner all the way to the output.
pub(crate) enum DiagramPair {
    Shared(FeynmanDiagram),
    Distinct {
        comparison: FeynmanDiagram,
        representative: Box<FeynmanDiagram>,
        reversed_edges: BTreeSet<EdgeId>,
    },
}

impl DiagramPair {
    pub(crate) fn comparison(&self) -> &FeynmanDiagram {
        match self {
            Self::Shared(diagram)
            | Self::Distinct {
                comparison: diagram,
                ..
            } => diagram,
        }
    }

    fn representative(&self) -> &FeynmanDiagram {
        match self {
            Self::Shared(diagram) => diagram,
            Self::Distinct { representative, .. } => representative,
        }
    }

    fn with_name(self, name: &str) -> Self {
        match self {
            Self::Shared(diagram) => Self::Shared(diagram.with_name(name)),
            Self::Distinct {
                comparison,
                representative,
                reversed_edges,
            } => Self::Distinct {
                comparison: comparison.with_name(name),
                representative: Box::new((*representative).with_name(name)),
                reversed_edges,
            },
        }
    }

    fn with_loop_momentum_edges(self, edges: &[EdgeId]) -> Result<Self, DiagramError> {
        match self {
            Self::Shared(diagram) => Ok(Self::Shared(diagram.with_loop_momentum_edges(edges)?)),
            Self::Distinct {
                comparison,
                representative,
                reversed_edges,
            } => Ok(Self::Distinct {
                comparison: comparison.with_loop_momentum_edges(edges)?,
                representative: Box::new((*representative).with_loop_momentum_edges(edges)?),
                reversed_edges,
            }),
        }
    }

    fn with_cuts(self, cuts: Vec<DiagramCut>) -> Result<Self, DiagramError> {
        match self {
            Self::Shared(diagram) => Ok(Self::Shared(diagram.with_cuts(cuts)?)),
            Self::Distinct {
                comparison,
                representative,
                reversed_edges,
            } => {
                let comparison = comparison.with_cuts(cuts)?;
                let partitions = transport_cut_metadata(comparison.cuts(), &reversed_edges)
                    .into_iter()
                    .map(|cut| (cut.left.half_edges, cut.right.half_edges))
                    .collect();
                Ok(Self::Distinct {
                    comparison,
                    representative: Box::new((*representative).with_cut_partitions(partitions)?),
                    reversed_edges,
                })
            }
        }
    }

    /// Read surviving physical source identities before merged partitions change their IDs.
    pub(crate) fn source_identity(
        &self,
        output: &GenerationOutput,
    ) -> Result<Option<SourceIdentity>, GenerationError> {
        if matches!(output, GenerationOutput::Minimal) {
            return Ok(None);
        }
        let source = self.representative();
        let name = source
            .name()
            .ok_or_else(|| GenerationError::InvalidGroups("named source has no name".into()))?;
        Ok(Some(SourceIdentity {
            id: source.id()?,
            name: name.to_owned(),
        }))
    }

    /// Grouping compares canonical charged flows. Move the physical owner to the output,
    /// transporting the comparison master's merged partitions, factor and routing onto it.
    pub(crate) fn into_representative(self) -> Result<FeynmanDiagram, DiagramError> {
        let (comparison, mut representative, reversed_edges) = match self {
            Self::Shared(diagram) => return Ok(diagram),
            Self::Distinct {
                comparison,
                representative,
                reversed_edges,
            } => (comparison, *representative, reversed_edges),
        };
        let cuts = transport_cut_metadata(comparison.cuts(), &reversed_edges);
        let thresholds = transport_threshold_candidate_metadata(
            comparison.topology_threshold_candidates(),
            &reversed_edges,
        );
        if cuts != representative.cuts() {
            representative = representative.with_cut_partitions(
                cuts.into_iter()
                    .map(|cut| (cut.left.half_edges, cut.right.half_edges))
                    .collect(),
            )?;
        }
        if thresholds != representative.topology_threshold_candidates() {
            representative = representative.with_topology_threshold_partitions(
                thresholds
                    .into_iter()
                    .map(|candidate| (candidate.left, candidate.right))
                    .collect(),
            )?;
        }
        if representative.loop_momentum_basis().loop_edges
            != comparison.loop_momentum_basis().loop_edges
        {
            representative = representative
                .with_loop_momentum_edges(&comparison.loop_momentum_basis().loop_edges)?;
        }
        if let Some(name) = comparison.name() {
            representative = representative.with_name(name);
        }
        Ok(representative.with_overall_factor(comparison.overall_factor().clone()))
    }
}

#[derive(Default)]
struct CutPartitionInventory {
    physical: Vec<DiagramCut>,
    topology_threshold_candidates: Vec<DiagramThresholdCandidate>,
}

fn transport_cut_half_edge(
    mut half_edge: DiagramHalfEdge,
    reversed_edges: &BTreeSet<EdgeId>,
) -> DiagramHalfEdge {
    if reversed_edges.contains(&half_edge.edge) {
        half_edge.endpoint = match half_edge.endpoint {
            DiagramEndpoint::Source => DiagramEndpoint::Target,
            DiagramEndpoint::Target => DiagramEndpoint::Source,
        };
    }
    half_edge
}

fn transport_cut_metadata(
    cuts: &[DiagramCut],
    reversed_edges: &BTreeSet<EdgeId>,
) -> Vec<DiagramCut> {
    cuts.iter()
        .map(|cut| DiagramCut {
            cut: cut
                .cut
                .iter()
                .copied()
                .map(|half_edge| transport_cut_half_edge(half_edge, reversed_edges))
                .collect(),
            left: DiagramCutSide {
                half_edges: cut
                    .left
                    .half_edges
                    .iter()
                    .copied()
                    .map(|half_edge| transport_cut_half_edge(half_edge, reversed_edges))
                    .collect(),
                coupling_orders: cut.left.coupling_orders.clone(),
                loop_count: cut.left.loop_count,
            },
            right: DiagramCutSide {
                half_edges: cut
                    .right
                    .half_edges
                    .iter()
                    .copied()
                    .map(|half_edge| transport_cut_half_edge(half_edge, reversed_edges))
                    .collect(),
                coupling_orders: cut.right.coupling_orders.clone(),
                loop_count: cut.right.loop_count,
            },
        })
        .collect()
}

fn transport_threshold_candidate_metadata(
    candidates: &[DiagramThresholdCandidate],
    reversed_edges: &BTreeSet<EdgeId>,
) -> Vec<DiagramThresholdCandidate> {
    candidates
        .iter()
        .map(|candidate| DiagramThresholdCandidate {
            cut: candidate
                .cut
                .iter()
                .copied()
                .map(|half_edge| transport_cut_half_edge(half_edge, reversed_edges))
                .collect(),
            left: candidate
                .left
                .iter()
                .copied()
                .map(|half_edge| transport_cut_half_edge(half_edge, reversed_edges))
                .collect(),
            right: candidate
                .right
                .iter()
                .copied()
                .map(|half_edge| transport_cut_half_edge(half_edge, reversed_edges))
                .collect(),
        })
        .collect()
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
struct FermionSigns {
    internal_loop_negative: bool,
    include_external_ordering: bool,
    external_ordering_negative: bool,
    include_antifermion_spin_sum: bool,
    antifermion_spin_sum_negative: bool,
}

struct NormalizedFermionFlow {
    graph: GeneratedGraph,
    reversed_edges: BTreeSet<EdgeId>,
    signs: FermionSigns,
}

struct FermionNormalizationState<'a> {
    visited_edges: &'a mut [bool],
    vertices: &'a mut [Option<(usize, usize)>],
    particles: &'a mut [Option<ParticleId>],
}

/// The fermion and ghost lines of a graph, in its plain-node numbering.
struct FermionLines<'a> {
    graph: &'a GeneratedGraph,
    nodes: PlainNodes,
    /// The plain nodes of the ends of every edge.
    ends: Vec<(usize, usize)>,
    /// For every node, the anticommuting edges at it with their other end.
    adjacency: Vec<Vec<(usize, usize)>>,
}

impl<'a> FermionLines<'a> {
    fn new(graph: &'a GeneratedGraph, model: &Model) -> Result<Self, GenerationError> {
        let nodes = graph.plain_nodes();
        let ends: Vec<_> = graph
            .edges()
            .enumerate()
            .map(|(id, edge)| nodes.ends(id, &edge))
            .collect();
        let mut adjacency = vec![Vec::new(); nodes.len()];
        let mut degrees = vec![0; nodes.len()];
        for (edge_id, edge) in graph.edges().enumerate() {
            let particle = model.particle_by_id(edge.particle)?;
            if !(particle.is_fermion() || particle.is_ghost()) {
                continue;
            }
            let (left, right) = ends[edge_id];
            adjacency[left].push((edge_id, right));
            if left == right {
                degrees[left] += 2;
            } else {
                adjacency[right].push((edge_id, left));
                degrees[left] += 1;
                degrees[right] += 1;
            }
        }
        if let Some((vertex, degree)) = degrees
            .iter()
            .copied()
            .enumerate()
            .find(|(_, degree)| *degree > 2)
        {
            return Err(GenerationError::UnsupportedFermionBranching { vertex, degree });
        }
        Ok(Self {
            graph,
            nodes,
            ends,
            adjacency,
        })
    }

    /// Orient the line through `starting_edge` and return the external legs it connects.
    fn normalize(
        &self,
        starting_edge: usize,
        state: &mut FermionNormalizationState<'_>,
        model: &Model,
        canonicalize_virtual_start: bool,
    ) -> Result<BTreeSet<usize>, GenerationError> {
        let (graph, nodes, ends, adjacency) =
            (self.graph, &self.nodes, &self.ends, &self.adjacency);
        let starting = graph.edge(starting_edge);
        let starting_particle = model.particle_by_id(starting.particle)?;
        let reverse_start = canonicalize_virtual_start && starting_particle.is_antiparticle();
        let (left, right) = if reverse_start {
            (ends[starting_edge].1, ends[starting_edge].0)
        } else {
            ends[starting_edge]
        };
        let chain_is_antiparticle = starting_particle.is_antiparticle() && !reverse_start;
        state.visited_edges[starting_edge] = true;
        state.vertices[starting_edge] = Some((left, right));
        state.particles[starting_edge] = Some(if reverse_start {
            starting_particle.antiparticle
        } else {
            starting.particle
        });
        let mut external_legs = BTreeSet::new();
        let mut record_external = |node: usize| {
            if let Some(external) = nodes.external(node) {
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
                .find(|(edge_id, _)| !state.visited_edges[*edge_id])
            {
                state.visited_edges[edge_id] = true;
                let vertices = if read_forward {
                    (current, next)
                } else {
                    (next, current)
                };
                state.vertices[edge_id] = Some(vertices);
                let edge = graph.edge(edge_id);
                let particle = model.particle_by_id(edge.particle)?;
                // Legacy GammaLoop represents a complete charged chain with the
                // same particle (or antiparticle) on every edge.  When a
                // canonical edge points against the chain, reversing only its
                // endpoints is not enough: choose the conjugate precisely when
                // needed to match the species selected by the starting edge.
                state.particles[edge_id] = Some(
                    if vertices != ends[edge_id]
                        && particle.is_antiparticle() != chain_is_antiparticle
                    {
                        particle.antiparticle
                    } else {
                        edge.particle
                    },
                );
                record_external(next);
                current = next;
            }
        }
        Ok(external_legs)
    }
}

fn all_incoming_particle<'model>(
    external: &ExternalNode,
    model: &'model Model,
) -> Result<&'model Particle, GenerationError> {
    let particle = model.particle_by_id(external.particle)?;
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

/// Intern exactly the reduced numerator colors in their old string order. External
/// particle names and internal mass/spin spellings share one table, including collisions.
/// A strictly order-preserving relabeling leaves the canonizer's partition and tie breaks
/// unchanged; the richer generator certificate has a different equivalence relation.
struct NumeratorColors {
    particles: Vec<usize>,
    masses_and_spins: Vec<usize>,
}

impl NumeratorColors {
    fn new(model: &Model) -> Result<Self, GenerationError> {
        let particles = model
            .particles()
            .iter()
            .map(|p| p.name.clone())
            .collect::<Vec<_>>();
        let masses_and_spins = model
            .particles()
            .iter()
            .map(|particle| {
                Ok(format!(
                    "{} | {}",
                    model.parameter_by_id(particle.mass)?.name,
                    particle.spin
                ))
            })
            .collect::<Result<Vec<_>, GenerationError>>()?;
        let ranks = particles
            .iter()
            .chain(&masses_and_spins)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .enumerate()
            .map(|(rank, text)| (text, rank))
            .collect::<BTreeMap<_, _>>();
        Ok(Self {
            particles: particles.iter().map(|text| ranks[text]).collect(),
            masses_and_spins: masses_and_spins.iter().map(|text| ranks[text]).collect(),
        })
    }

    fn edge(
        &self,
        edge: &GeneratedEdge,
        model: &Model,
        exact: bool,
    ) -> Result<usize, GenerationError> {
        let particle = model.particle_by_id(edge.particle)?;
        let id = if edge.is_external() && edge.sink.external().is_some() {
            particle.antiparticle
        } else {
            edge.particle
        };
        Ok(if edge.is_external() || exact {
            self.particles[id.index()]
        } else {
            self.masses_and_spins[id.index()]
        })
    }
}

struct ResolvedProcess {
    numerator_colors: NumeratorColors,
    generation_type: GenerationType,
    incoming: Vec<Particle>,
    outgoing: Vec<Vec<Particle>>,
    symmetrize_initial: bool,
    symmetrize_final: bool,
    symmetrize_left_right: bool,
    symmetrize_external_fermions: bool,
}

impl ResolvedProcess {
    fn new(
        model: &Model,
        process: &Process,
        generation_type: GenerationType,
        options: &GenerationOptions,
    ) -> Result<Self, GenerationError> {
        let resolve = |selector: &ParticleSelector| -> Result<Particle, GenerationError> {
            Ok(model.particle_by_id(selector.resolve(model)?)?.clone())
        };
        Ok(Self {
            numerator_colors: NumeratorColors::new(model)?,
            generation_type,
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
            symmetrize_initial: options.symmetrizes_initial(),
            symmetrize_final: options.symmetrizes_final(),
            symmetrize_left_right: options.symmetrizes_left_right(),
            symmetrize_external_fermions: options.symmetrizes_external_fermions(),
        })
    }

    /// The external legs of the generated graphs: incoming legs first, then the outgoing legs
    /// (the incoming states again for a cross section).
    ///
    /// The outgoing legs of a cross section are the sewn copies of the incoming states, the
    /// initial state of the conjugate amplitude, so only `symmetrize_initial` applies to them:
    /// exchanging two identical beams must move each sewn copy with its incoming leg. Such a
    /// diagonal permutation is not a symmetry of independent leg labels, so every cross-section
    /// leg is labelled exactly and the interchangeable sewn pairs are identified after
    /// generation; see [`ProcessExternals::sewn_pair_orbit`]. The exchange is a rotation of the
    /// beams only in a `2 -> X` process, so identical beams among more incoming particles are
    /// an error.
    fn externals(&self, model: &Model) -> Result<ProcessExternals, GenerationError> {
        let incoming = self
            .incoming
            .iter()
            .map(|particle| (particle, ExternalState::Incoming));
        let outgoing = match self.generation_type {
            GenerationType::Amplitude => &self.outgoing[0],
            GenerationType::CrossSection => &self.incoming,
        }
        .iter()
        .map(|particle| (particle, ExternalState::Outgoing));
        let externals = incoming
            .chain(outgoing)
            .enumerate()
            .map(|(index, (particle, state))| self.external(model, particle, index, state))
            .collect::<Result<Vec<_>, _>>()?;
        let mut sewn_pairs = BTreeMap::<&str, Vec<_>>::new();
        if self.generation_type == GenerationType::CrossSection && self.symmetrize_initial {
            for (index, particle) in self.incoming.iter().enumerate() {
                sewn_pairs.entry(&particle.name).or_default().push([
                    ExternalLabel::Exact {
                        state: ExternalState::Incoming,
                        index,
                    },
                    ExternalLabel::Exact {
                        state: ExternalState::Outgoing,
                        index: self.incoming.len() + index,
                    },
                ]);
            }
        }
        let sewn_pairs = sewn_pairs
            .into_values()
            .filter(|pairs| pairs.len() > 1)
            .collect::<Vec<_>>();
        if !sewn_pairs.is_empty() && self.incoming.len() != 2 {
            return Err(GenerationError::UnsupportedBeamSymmetrization {
                incoming: self.incoming.len(),
            });
        }
        Ok(ProcessExternals::new(externals, sewn_pairs))
    }

    fn external(
        &self,
        model: &Model,
        particle: &Particle,
        index: usize,
        state: ExternalState,
    ) -> Result<(ExternalLabel, ParticleStub, ExternalNode), GenerationError> {
        // An outgoing particle enters the graph as its antiparticle.
        let entering = match state {
            ExternalState::Incoming => particle,
            ExternalState::Outgoing => model.antiparticle(particle)?,
        };
        let symmetrized_within_state = self.generation_type == GenerationType::Amplitude
            && match state {
                ExternalState::Incoming => self.symmetrize_initial,
                ExternalState::Outgoing => self.symmetrize_final,
            }
            && (!particle.is_fermion() || self.symmetrize_external_fermions);
        // Left/right symmetry is one global transformation of a sewn graph,
        // not an independent automorphism of every external pair. It is
        // applied after generation when canonical diagram keys are built.
        let label = if symmetrized_within_state {
            ExternalLabel::WithinState(state)
        } else {
            ExternalLabel::Exact { state, index }
        };
        Ok((
            label,
            incoming_stub(model, entering),
            ExternalNode {
                index,
                state,
                particle: model.particle_id(&particle.name)?,
            },
        ))
    }

    /// Whether a sewn graph and its mirror image, which exchanges the amplitude sides, are
    /// identified.
    fn identifies_mirror_images(&self) -> bool {
        self.generation_type == GenerationType::CrossSection && self.symmetrize_left_right
    }

    /// Choose the reduced-skeleton representative before any numerator
    /// indices are instantiated.
    ///
    /// For left/right symmetry the extra orbit operation is global: either
    /// every sewn incoming leg is exchanged with its outgoing partner, or none
    /// is. Keeping the exact external colors in each skeleton prevents graph
    /// automorphisms from swapping only a subset of the pairs. Once the
    /// smallest skeleton is selected, the colored graph adopts
    /// that canonical vertex/edge frame so UFO slot matching and propagator
    /// momenta are instantiated identically for isomorphic diagrams.
    fn canonicalize_numerator_graph(
        &self,
        graph: &GeneratedGraph,
        model: &Model,
        grouping: &NumeratorGrouping,
    ) -> Result<GeneratedGraph, GenerationError> {
        let grouping_options = match grouping {
            NumeratorGrouping::Identical(options)
            | NumeratorGrouping::UpToSign(options)
            | NumeratorGrouping::UpToScalar(options) => Some(options),
            NumeratorGrouping::None => None,
        };
        let include_left_right = self.identifies_mirror_images();
        let exact_internal_particles =
            grouping_options.is_some_and(|options| !options.differentiate_particle_masses_only);

        let selected = self.canonicalize_numerator_graph_pass(
            graph,
            model,
            exact_internal_particles,
            include_left_right,
            include_left_right,
        )?;
        let canonical = if include_left_right {
            self.canonicalize_numerator_graph_pass(
                &selected,
                model,
                exact_internal_particles,
                false,
                true,
            )?
        } else {
            selected
        };
        Ok(canonical)
    }

    fn canonicalize_numerator_graph_pass(
        &self,
        graph: &GeneratedGraph,
        model: &Model,
        exact_internal_particles: bool,
        include_left_right_remapping: bool,
        use_left_right_external_colors: bool,
    ) -> Result<GeneratedGraph, GenerationError> {
        let incoming_count = self.incoming.len();
        let externals: BTreeMap<usize, ExternalNode> = if include_left_right_remapping {
            graph
                .externals()
                .map(|(_, external, _)| (external.index, external))
                .collect()
        } else {
            BTreeMap::new()
        };
        if include_left_right_remapping && externals.len() != 2 * incoming_count {
            return Err(GenerationError::GraphConstruction(format!(
                "left/right canonicalization expected {} sewn external legs, found {}",
                2 * incoming_count,
                externals.len()
            )));
        }
        // The global left/right exchange swaps every sewn incoming leg with its outgoing
        // partner.
        let partner = |external: ExternalNode| {
            let index = if external.index < incoming_count {
                external.index + incoming_count
            } else {
                external.index - incoming_count
            };
            externals[&index]
        };
        let external_color = |external: ExternalNode| -> Result<i64, GenerationError> {
            let overflow = || GenerationError::ArithmeticOverflow("an external canonical color");
            let index = i64::try_from(external.index).map_err(|_| overflow())?;
            let incoming_count = i64::try_from(incoming_count).map_err(|_| overflow())?;
            Ok(if use_left_right_external_colors {
                if index < incoming_count {
                    -2001 - index
                } else {
                    -3001 - (index - incoming_count)
                }
            } else {
                index + 1
            })
        };
        // External edges are colored by the particle entering the interaction graph.
        // Internal direction is omitted from the reduced skeleton but restored below.
        // The two global left/right choices share incidence and edge colors. Prepare that
        // skeleton once and vary only its external half-edge keys when canonizing it.
        let mut skeleton = HedgeGraphBuilder::new();
        let nodes: Vec<NodeIndex> = (0..graph.n_vertices())
            .map(|_| skeleton.add_node(()))
            .collect();
        for edge in graph.edges() {
            let color = self
                .numerator_colors
                .edge(&edge, model, exact_internal_particles)?;
            match (edge.source, edge.sink) {
                (Endpoint::Vertex(source), Endpoint::Vertex(sink)) => {
                    skeleton.add_edge(nodes[source], nodes[sink], color, Orientation::Undirected)
                }
                (Endpoint::External(external), Endpoint::Vertex(vertex))
                | (Endpoint::Vertex(vertex), Endpoint::External(external)) => {
                    let normal = external_color(external)?;
                    let swapped = if include_left_right_remapping {
                        external_color(partner(external))?
                    } else {
                        normal
                    };
                    skeleton.add_external_edge(
                        nodes[vertex].add_data(Some([normal, swapped])),
                        color,
                        Orientation::Default,
                        Flow::Sink,
                    )
                }
                (Endpoint::External(_), Endpoint::External(_)) => {
                    unreachable!("generated edges touch an interaction vertex")
                }
            }
        }
        let skeleton: HedgeGraph<usize, (), Option<[i64; 2]>> = skeleton.build();
        let mut candidates = Vec::new();
        for swapped in [false, true] {
            if swapped && !include_left_right_remapping {
                continue;
            }
            let canonical = skeleton.canonize_by(
                ignore,
                |color| *color,
                |colors| colors.map(|colors| colors[usize::from(swapped)]),
            );
            let positions: Vec<usize> = (0..graph.n_vertices())
                .map(|vertex| canonical.canonical_node(NodeIndex(vertex)).0)
                .collect();
            candidates.push((canonical.into_canonical_form(), swapped, positions));
        }
        let (_, left_right_swapped, positions) = candidates
            .into_iter()
            .min_by(|left, right| (&left.0, left.1).cmp(&(&right.0, right.1)))
            .ok_or_else(|| GenerationError::GraphConstruction("no canonical graph".to_owned()))?;

        // The selected global exchange is an involution. Applying it to the legs keeps each
        // physical external label attached to its original process state while adopting the
        // selected topology.
        let endpoint = |endpoint: Endpoint| match endpoint {
            Endpoint::Vertex(vertex) => Endpoint::Vertex(positions[vertex]),
            Endpoint::External(external) if left_right_swapped => {
                Endpoint::External(partner(external))
            }
            Endpoint::External(external) => Endpoint::External(external),
        };
        let mut edges = graph
            .edges()
            .enumerate()
            .map(|(edge_id, edge)| {
                // All external edges point into the interaction graph and internal
                // edges point from the lower to the higher canonical vertex.
                let flipped = edge.sink.external().is_some()
                    || matches!(
                        (edge.source, edge.sink),
                        (Endpoint::Vertex(source), Endpoint::Vertex(sink))
                            if positions[source] > positions[sink]
                    );
                let mut particle = edge.particle;
                if edge.directed && flipped != left_right_swapped {
                    particle = model.particle_by_id(particle)?.antiparticle;
                }
                let (source, sink) = if flipped {
                    (edge.sink, edge.source)
                } else {
                    (edge.source, edge.sink)
                };
                let pdg_code = model.particle_by_id(particle)?.pdg_code;
                Ok((
                    GeneratedEdge {
                        source: endpoint(source),
                        sink: endpoint(sink),
                        particle,
                        directed: edge.directed,
                    },
                    pdg_code,
                    edge_id,
                    flipped,
                ))
            })
            .collect::<Result<Vec<_>, GenerationError>>()?;
        // Match the UFO/GammaLoop canonical momentum convention, which diagram finalization
        // also uses for internal edges. Signed PDG order is deliberately used here rather than
        // model insertion order: CP-conjugate charged lines can otherwise acquire the same
        // edge/LMB order even when their legacy canonical representatives differ.
        edges.sort_by_key(|(edge, pdg_code, _, _)| (edge.source, edge.sink, *pdg_code));
        Ok(graph.reframed(
            &positions,
            edges
                .into_iter()
                .map(|(edge, _, edge_id, flipped)| (edge_id, flipped, edge.particle)),
            |external| {
                if left_right_swapped {
                    partner(external)
                } else {
                    external
                }
            },
        ))
    }

    fn repair_cp_interaction_rules(
        graph: GeneratedGraph,
        model: &Model,
    ) -> Result<GeneratedGraph, GenerationError> {
        let mut replacements = BTreeMap::new();
        for (vertex, rule_id) in graph.vertices().enumerate() {
            let rule = model.vertex_rule_by_id(rule_id)?;
            let expected = vertex_rule_signature(model, rule)?;
            let mut actual = available_interaction_half_edges(&graph, vertex, model)?
                .into_iter()
                .map(|(_, signature)| signature)
                .collect::<Vec<_>>();
            actual.sort();
            if actual == expected {
                continue;
            }

            let mut cp_expected = expected
                .iter()
                .cloned()
                .map(ParticleStub::flipped)
                .collect::<Vec<_>>();
            cp_expected.sort();
            let coupling_orders = rule.coupling_orders(model);
            let mut candidates = Vec::new();
            if actual == cp_expected {
                for (candidate_index, candidate) in model.vertex_rules().iter().enumerate() {
                    if candidate.coupling_orders(model) == coupling_orders
                        && vertex_rule_signature(model, candidate)? == actual
                    {
                        candidates.push((
                            model.vertex_rule_id_at(candidate_index)?,
                            candidate.name.clone(),
                        ));
                    }
                }
            }
            if candidates.len() != 1 {
                return Err(GenerationError::InvalidCpInteractionRule {
                    vertex,
                    interaction: rule.name.clone(),
                    candidates: candidates.into_iter().map(|(_, name)| name).collect(),
                });
            }
            replacements.insert(vertex, candidates[0].0);
        }

        if replacements.is_empty() {
            return Ok(graph);
        }
        let rules = graph
            .vertices()
            .enumerate()
            .map(|(vertex, rule)| replacements.get(&vertex).copied().unwrap_or(rule))
            .collect::<Vec<_>>();
        Ok(graph.with_vertex_rules(rules))
    }

    // Normalize every fermion chain before exposing the diagram. A chain keeps
    // the orientation of its first external edge and all following edges are
    // oriented consistently with it. Purely virtual loops start from their
    // first edge in canonical edge order.
    fn normalize_fermion_flows(
        &self,
        graph: &GeneratedGraph,
        model: &Model,
    ) -> Result<NormalizedFermionFlow, GenerationError> {
        let lines = FermionLines::new(graph, model)?;
        let (nodes, ends) = (&lines.nodes, &lines.ends);
        let external_nodes: BTreeMap<_, _> = (0..nodes.len())
            .filter_map(|node| {
                nodes
                    .external(node)
                    .map(|external| (external.index, (node, external)))
            })
            .collect();
        let mut visited_edges = vec![false; graph.n_edges()];
        let mut normalized_vertices = vec![None; graph.n_edges()];
        let mut normalized_particles = vec![None; graph.n_edges()];
        let mut normalization = FermionNormalizationState {
            visited_edges: &mut visited_edges,
            vertices: &mut normalized_vertices,
            particles: &mut normalized_particles,
        };
        // Pairing of the external fermion flows. The keys are sorted pairs of
        // all-incoming PDGs, and values are connected external leg IDs in the
        // order of that key.
        let mut pairings: BTreeMap<(i64, i64), Vec<(usize, usize)>> = BTreeMap::new();
        let mut line_pairings = Vec::new();

        // First fix flows connected to external legs, then all internal flows.
        for (node_id, external) in external_nodes.values() {
            let Some((edge_id, _)) = lines.adjacency[*node_id]
                .iter()
                .find(|(edge_id, _)| !normalization.visited_edges[*edge_id])
            else {
                continue;
            };
            let legs = lines.normalize(*edge_id, &mut normalization, model, false)?;
            let legs = legs.into_iter().collect::<Vec<_>>();
            // Both spinor and ghost chains anticommute. Ghosts therefore
            // enter external permutation signs even though their scalar
            // propagators do not contribute antifermion spin sums.
            if !model
                .particle_by_id(external.particle)?
                .is_anticommutating()
            {
                continue;
            }
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
                first_particle.is_anticommutating(),
                first_particle.is_antiparticle(),
                second_particle.is_anticommutating(),
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
        for edge_id in 0..graph.n_edges() {
            let particle = model.particle_by_id(graph.edge(edge_id).particle)?;
            if (particle.is_fermion() || particle.is_ghost())
                && !normalization.visited_edges[edge_id]
            {
                lines.normalize(edge_id, &mut normalization, model, true)?;
            }
        }

        let mut reversed_edges = BTreeSet::new();
        let mut edges = Vec::with_capacity(graph.n_edges());
        for (edge_id, edge) in graph.edges().enumerate() {
            let particle = model.particle_by_id(edge.particle)?;
            let canonicalize_virtual_antiparticle = normalized_vertices[edge_id].is_none()
                && !edge.is_external()
                && particle.is_antiparticle()
                && !(particle.is_fermion() || particle.is_ghost());
            let (left, right) = normalized_vertices[edge_id].unwrap_or({
                if canonicalize_virtual_antiparticle {
                    (ends[edge_id].1, ends[edge_id].0)
                } else {
                    ends[edge_id]
                }
            });
            let particle = normalized_particles[edge_id].unwrap_or({
                if canonicalize_virtual_antiparticle {
                    particle.antiparticle
                } else {
                    edge.particle
                }
            });
            let reversed = (left, right) != ends[edge_id];
            if reversed {
                reversed_edges.insert(EdgeId(edge_id));
            }
            edges.push((edge_id, reversed, particle));
        }
        let positions = (0..graph.n_vertices()).collect::<Vec<_>>();
        let normalized_graph = graph.reframed(&positions, edges, |external| external);

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
                // Closing n open fermion chains into c cycles has permutation
                // parity (-1)^(n-c). The closure itself is not a virtual loop.
                (line_pairings.len() + self.external_fermion_loop_count(&line_pairings)?) % 2 == 1
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
            reversed_edges,
            signs: FermionSigns {
                internal_loop_negative: false,
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
            if particle.is_anticommutating() {
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

    // Enumerate exact s-channel cuts for cross-section generation. The caller
    // applies amplitude-level constraints to the returned partitions.
    fn cut_partitions(
        &self,
        graph: &GeneratedGraph,
        model: &Model,
        options: &GenerationOptions,
    ) -> Result<CutPartitionInventory, GenerationError> {
        let unresolved = unresolved_cut_content(model, options)?;
        let targets: Vec<Vec<ParticleId>> = self
            .outgoing
            .iter()
            .map(|particles| {
                let mut particle_ids = particles
                    .iter()
                    .map(|particle| model.particle_id(&particle.name))
                    .collect::<Result<Vec<_>, _>>()?;
                particle_ids.sort();
                Ok(particle_ids)
            })
            .collect::<Result<_, ModelError>>()?;
        let he_graph = graph.with_external_vertices();
        // Legacy cross-section cuts were enumerated after sewing each incoming
        // external leg to its outgoing partner.  The resulting source/target
        // nodes are the internal attachment nodes, and an initial-state edge
        // can therefore never itself be part of an s-channel cut.  On the
        // canonical pre-sewn FeynKit graph, force every degree-one external
        // vertex and its internal attachment onto the same side to reproduce
        // that topology exactly.
        let mut source = BTreeSet::new();
        let mut target = BTreeSet::new();
        for (node_id, neighbours, node) in he_graph.iter_nodes() {
            if let PlainVertex::External(external) = node {
                let incident = neighbours.collect::<Vec<_>>();
                if incident.len() != 1 {
                    return Err(GenerationError::GraphConstruction(format!(
                        "cross-section external vertex {node_id:?} has degree {}, expected one",
                        incident.len(),
                    )));
                }
                let attachment = he_graph.node_id(he_graph.inv(incident[0]));
                match external.state {
                    ExternalState::Incoming => {
                        source.insert(node_id);
                        source.insert(attachment);
                    }
                    ExternalState::Outgoing => {
                        target.insert(node_id);
                        target.insert(attachment);
                    }
                }
            }
        }
        if source.is_empty() || target.is_empty() {
            // Vacuum or degenerate graphs have no physical s-channel cut.
            return Ok(CutPartitionInventory::default());
        }
        let initial_state_edges = graph
            .edges()
            .enumerate()
            .filter_map(|(edge, data)| data.is_external().then_some(EdgeId(edge)))
            .collect::<BTreeSet<_>>();

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

        let stable_half_edges = |subgraph: &SuBitGraph| {
            let mut half_edges = subgraph
                .included_iter()
                .map(|hedge| DiagramHalfEdge {
                    edge: EdgeId(he_graph[&hedge].0),
                    endpoint: match he_graph.flow(hedge) {
                        Flow::Source => DiagramEndpoint::Source,
                        Flow::Sink => DiagramEndpoint::Target,
                    },
                })
                .collect::<Vec<_>>();
            half_edges.sort();
            half_edges.dedup();
            half_edges
        };
        let summarize_side = |subgraph: &SuBitGraph| -> Result<DiagramCutSide, GenerationError> {
            let mut coupling_orders = BTreeMap::new();
            for (_, _, node) in he_graph.iter_nodes_of(subgraph) {
                let PlainVertex::Interaction(rule_id) = node else {
                    continue;
                };
                for (name, order) in model.vertex_rule_by_id(*rule_id)?.coupling_orders(model) {
                    *coupling_orders.entry(name).or_insert(0) += order;
                }
            }
            let internal = InternalSubGraph::cleaned_filter_pessimist(subgraph.clone(), &he_graph);
            Ok(DiagramCutSide {
                half_edges: stable_half_edges(subgraph),
                coupling_orders,
                loop_count: he_graph.cyclotomatic_number(&internal),
            })
        };

        let mut inventory = CutPartitionInventory::default();
        let source = source.into_iter().collect::<Vec<_>>();
        let target = target.into_iter().collect::<Vec<_>>();
        for (left, cut, right) in he_graph.all_cuts_from_ids(&source, &target) {
            let cut_half_edges = stable_half_edges(&cut.left);
            if cut_half_edges
                .iter()
                .any(|half_edge| initial_state_edges.contains(&half_edge.edge))
            {
                // Keep this guard even though the seeded attachment nodes make
                // such a cut impossible: serialized threshold candidates must
                // never rely on the runtime bridge stripping a crossing edge.
                continue;
            }
            if cut_half_edges
                .iter()
                .filter(|half_edge| !initial_state_edges.contains(&half_edge.edge))
                .count()
                > 1
            {
                inventory
                    .topology_threshold_candidates
                    .push(DiagramThresholdCandidate {
                        cut: cut_half_edges.clone(),
                        left: stable_half_edges(&left),
                        right: stable_half_edges(&right),
                    });
            }
            if !validate_connectivity(&left) || !validate_connectivity(&right) {
                continue;
            }
            // Convert sink-oriented cut half-edges to antiparticles so the
            // content describes the physical state crossing left to right.
            let mut cut_content = cut
                .iter_left_hedges()
                .map(|hedge| {
                    let particle = he_graph[[&hedge]];
                    if matches!(he_graph.flow(hedge), Flow::Sink) {
                        Ok(model.particle_by_id(particle)?.antiparticle)
                    } else {
                        Ok(particle)
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
            inventory.physical.push(DiagramCut {
                cut: cut_half_edges,
                left: summarize_side(&left)?,
                right: summarize_side(&right)?,
            });
        }
        // The cuts arrive in no fixed order, which finalization fixes for every diagram. Merging
        // mirror images compares inventories before that, so equal graphs need equal ones.
        if self.identifies_mirror_images() {
            inventory.physical.sort();
            inventory.topology_threshold_candidates.sort();
        }
        Ok(inventory)
    }

    fn matches_cut_content(
        cut_content: &[ParticleId],
        target: &[ParticleId],
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

    fn passes_sewn_filter(&self, graph: &GeneratedGraph, options: &GenerationOptions) -> bool {
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

        let mut graph = graph.with_external_vertices();
        let external_nodes: Vec<_> = graph
            .iter_nodes()
            .filter_map(|(node, _, color)| match color {
                PlainVertex::External(external) => Some((external.index, node, *color)),
                PlainVertex::Interaction(_) => None,
            })
            .collect();
        let externals: Vec<_> = external_nodes.iter().map(|(_, node, _)| *node).collect();
        if externals.is_empty() {
            return true;
        }
        let connected_components_before = graph
            .tadpoles(&externals)
            .expect("external nodes are distinct graph terminals with incident edges")
            .len()
            + 1;
        let by_index: BTreeMap<usize, (NodeIndex, PlainVertex)> = external_nodes
            .into_iter()
            .map(|(index, node, color)| (index, (node, color)))
            .collect();
        for incoming in 0..self.incoming.len() {
            let outgoing = self.incoming.len() + incoming;
            if let (Some((left, color)), Some((right, _))) =
                (by_index.get(&incoming), by_index.get(&outgoing))
            {
                graph.identify_nodes(&[*left, *right], *color);
            }
        }
        // Identification keeps the absorbed nodes until compaction; cycle
        // detection runs on the compacted quotient graph.
        graph.compact();
        let non_bridges = graph.non_bridges();
        graph.count_connected_components(&non_bridges) == connected_components_before
    }
}

#[cfg(test)]
mod tests {
    use crate::{CancellationToken, GraphGroupingOptions, NumeratorGrouping, VertexSelector};
    use feynkit_model::Model;
    use symbolica::parser::ParseSettings;

    use super::*;

    fn test_atom(expression: &str) -> Atom {
        Atom::parse(expression, "UFO", ParseSettings::default()).unwrap()
    }

    fn factor_atom(expression: &str) -> Atom {
        Atom::parse(
            expression,
            "feynkit_generator_factor",
            ParseSettings::default(),
        )
        .unwrap()
    }

    fn function_count(expression: &Atom, head: symbolica::atom::Symbol) -> usize {
        let mut count = 0;
        let _ = expression.replace_map(|term, _, _| {
            if matches!(term, AtomView::Fun(function) if function.get_symbol() == head) {
                count += 1;
            }
        });
        count
    }

    fn tensor_test_legs() -> [NumeratorHalfEdge; 4] {
        [
            NumeratorHalfEdge {
                edge: 10,
                flow: Flow::Source,
                hedge: 20,
                spin: 2,
                color: -3,
            },
            NumeratorHalfEdge {
                edge: 20,
                flow: Flow::Sink,
                hedge: 41,
                spin: 2,
                color: 3,
            },
            NumeratorHalfEdge {
                edge: 30,
                flow: Flow::Source,
                hedge: 60,
                spin: 3,
                color: 8,
            },
            NumeratorHalfEdge {
                edge: 40,
                flow: Flow::Sink,
                hedge: 81,
                spin: 3,
                color: 8,
            },
        ]
    }

    fn scalar_model() -> Model {
        Model::from_json(
            r#"{
                "name":"phi3",
                "restriction":null,
                "orders":[{"name":"QED","expansion_order":99,"hierarchy":1}],
                "parameters":[
                    {"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null},
                    {"name":"M","lhablock":"MASS","lhacode":[25],"nature":"external","parameter_type":"real","value":[1.0,0.0],"expression":null}
                ],
                "particles":[
                    {"pdg_code":25,"name":"phi","antiname":"phi","spin":1,"color":1,"mass":"M","width":"ZERO","texname":"phi","antitexname":"phi","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0}
                ],
                "propagators":[
                    {"name":"phi_prop","particle":"phi","numerator":"1","denominator":"P^2-M^2"}
                ],
                "lorentz_structures":[
                    {"name":"L1","spins":[1,1,1],"structure":"1"}
                ],
                "couplings":[
                    {"name":"GC1","expression":"g","orders":[["QED",1]],"value":null}
                ],
                "vertex_rules":[
                    {"name":"V1","particles":["phi","phi","phi"],"color_structures":["1"],"lorentz_structures":["L1"],"couplings":[["GC1"]]}
                ],
                "functions":[],
                "form_factors":[]
            }"#,
        )
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
                "parameters":[
                    {"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal",
                     "parameter_type":"real","value":[0.0,0.0],"expression":null},
                    {"name":"M","lhablock":"MASS","lhacode":[25],"nature":"external",
                     "parameter_type":"real","value":[1.0,0.0],"expression":null}
                ],
                "particles":[
                    {"pdg_code":1,"name":"f","antiname":"f~","spin":2,"color":1,
                     "mass":"ZERO","width":"ZERO","texname":"f","antitexname":"f~",
                     "charge":-1.0,"ghost_number":0,"lepton_number":1,"y_charge":0},
                    {"pdg_code":-1,"name":"f~","antiname":"f","spin":2,"color":1,
                     "mass":"ZERO","width":"ZERO","texname":"f~","antitexname":"f",
                     "charge":1.0,"ghost_number":0,"lepton_number":-1,"y_charge":0},
                    {"pdg_code":22,"name":"a","antiname":"a","spin":3,"color":1,
                     "mass":"ZERO","width":"ZERO","texname":"a","antitexname":"a",
                     "charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":25,"name":"phi","antiname":"phi","spin":1,"color":1,
                     "mass":"M","width":"ZERO","texname":"phi","antitexname":"phi",
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
                    "color_structures":["1","T(3,2,1)"],
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
                     "charge":"-1/3","ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":-5,"name":"b~","antiname":"b","spin":2,"color":-3,
                     "mass":"M","width":"ZERO","texname":"b~","antitexname":"b",
                     "charge":"1/3","ghost_number":0,"lepton_number":0,"y_charge":0},
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

    fn standard_model() -> Model {
        Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap()
    }

    fn photon_process() -> Process {
        Process::new([22_i64], [22_i64])
    }

    fn scalar_cross_section_process() -> Process {
        Process::new(["phi", "phi"], ["phi", "phi"])
    }

    fn particle(model: &Model, pdg: i64) -> ParticleId {
        model.particle_id_by_pdg(pdg).unwrap()
    }

    fn vertex(model: &Model, name: &str) -> VertexRuleId {
        model.vertex_rule_id(name).unwrap()
    }

    fn edge(model: &Model, pdg: i64) -> ParticleId {
        particle(model, pdg)
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum ColoredNode {
        External(ExternalNode),
        Interaction(VertexRuleId),
    }

    /// A generated-graph fixture written with degree-one external vertices; [`Self::build`]
    /// turns them into legs and keeps the edge order.
    #[derive(Debug, Clone, Default)]
    struct Graph {
        nodes: Vec<ColoredNode>,
        edges: Vec<(usize, usize, bool, ParticleId)>,
    }

    impl Graph {
        fn new() -> Self {
            Self::default()
        }

        fn add_node(&mut self, node: ColoredNode) -> usize {
            self.nodes.push(node);
            self.nodes.len() - 1
        }

        /// Undirected edges are stored from the lower to the higher node, as Graphica did.
        fn add_edge(
            &mut self,
            source: usize,
            sink: usize,
            directed: bool,
            particle: ParticleId,
        ) -> Result<usize, String> {
            let (source, sink) = if !directed && source > sink {
                (sink, source)
            } else {
                (source, sink)
            };
            self.edges.push((source, sink, directed, particle));
            Ok(self.edges.len() - 1)
        }

        fn interactions(&self) -> Vec<usize> {
            (0..self.nodes.len())
                .filter(|node| matches!(self.nodes[*node], ColoredNode::Interaction(_)))
                .collect()
        }

        fn build(&self) -> GeneratedGraph {
            let interactions = self.interactions();
            let endpoint = |node: usize| match self.nodes[node] {
                ColoredNode::External(external) => Endpoint::External(external),
                ColoredNode::Interaction(_) => Endpoint::Vertex(
                    interactions
                        .iter()
                        .position(|vertex| *vertex == node)
                        .unwrap(),
                ),
            };
            GeneratedGraph::new(
                interactions.iter().map(|node| match self.nodes[*node] {
                    ColoredNode::Interaction(rule) => rule,
                    ColoredNode::External(_) => unreachable!(),
                }),
                self.edges
                    .iter()
                    .map(|&(source, sink, directed, particle)| GeneratedEdge {
                        source: endpoint(source),
                        sink: endpoint(sink),
                        particle,
                        directed,
                    }),
            )
        }

        /// The fixture nodes of the ends of `edge` in a graph built from this fixture.
        fn ends(&self, graph: &GeneratedGraph, edge: usize) -> (usize, usize) {
            let interactions = self.interactions();
            let node = |endpoint: Endpoint| match endpoint {
                Endpoint::Vertex(vertex) => interactions[vertex],
                Endpoint::External(external) => self
                    .nodes
                    .iter()
                    .position(|node| *node == ColoredNode::External(external))
                    .unwrap(),
            };
            let edge = graph.edge(edge);
            (node(edge.source), node(edge.sink))
        }
    }

    #[test]
    fn loop_one_particle_irreducibility_preserves_trees_in_mixed_ranges() {
        let generator = Generator::new(scalar_model());
        let process = Process::new(["phi", "phi"], ["phi", "phi"]);
        let options = GenerationOptions::default()
            .with_loop_count(0, 1)
            .unwrap()
            .threads(1)
            .max_vertices(4);
        let automatic = options
            .clone()
            .with_graph_filter(GenerationFilter::LoopOneParticleIrreducible);
        assert_eq!(automatic.max_bridges(0), None);
        assert_eq!(automatic.max_bridges(1), Some(0));

        let all = generator
            .generate(&process, &options, GenerationType::Amplitude)
            .unwrap();
        let filtered = generator
            .generate(&process, &automatic, GenerationType::Amplitude)
            .unwrap();
        let strict = generator
            .generate(
                &process,
                &options.with_graph_filter(GenerationFilter::MaxNumberOfBridges(0)),
                GenerationType::Amplitude,
            )
            .unwrap();
        let ids = |result: &GenerationResult, loops| {
            result
                .diagrams
                .iter()
                .filter(|diagram| diagram.loop_count() == loops)
                .map(|diagram| diagram.id().unwrap())
                .collect::<BTreeSet<_>>()
        };
        assert_eq!(
            ids(&filtered, 0).len(),
            3,
            "all tree exchange channels survive"
        );
        assert_eq!(ids(&filtered, 0), ids(&all, 0));
        assert!(
            ids(&strict, 0).is_empty(),
            "an explicit bridge limit remains literal"
        );
        assert!(!ids(&filtered, 1).is_empty());
        assert_eq!(ids(&filtered, 1), ids(&strict, 1));
        assert!(
            ids(&all, 1).len() > ids(&filtered, 1).len(),
            "1PR loop diagrams are removed"
        );

        let loops_only = generator
            .generate(
                &process,
                &automatic.with_loop_count(1, 1).unwrap(),
                GenerationType::Amplitude,
            )
            .unwrap();
        assert_eq!(
            ids(&loops_only, 1),
            ids(&filtered, 1),
            "early pruning preserves the result"
        );
    }

    #[test]
    fn generates_tree_and_one_loop_diagrams() {
        let generator = Generator::new(scalar_model());
        let process = Process::new(["phi"], ["phi", "phi"]);
        let generated = generator
            .generate(
                &process,
                &GenerationOptions::default()
                    .with_loop_count(0, 1)
                    .unwrap()
                    .allow_self_loops(true),
                GenerationType::Amplitude,
            )
            .unwrap();
        assert!(!generated.diagrams.is_empty());
        assert!(generated.report.completed);
        assert!(
            generated
                .diagrams
                .iter()
                .all(|diagram| !diagram.numerator().is_zero())
        );
    }

    #[test]
    fn finalizes_external_fermion_projectors_unless_explicitly_overridden() {
        let generator = Generator::new(fermion_model());
        let process = Process::new([1_i64, -1], [1_i64, -1]);
        let options = GenerationOptions::default()
            .with_loop_count(0, 0)
            .unwrap()
            .threads(1)
            .max_vertices(2);
        let generated = generator
            .generate(&process, &options, GenerationType::Amplitude)
            .unwrap();
        assert!(!generated.diagrams.is_empty());
        for diagram in &generated.diagrams {
            for wavefunction in [
                feynkit_graph::symbols::u(),
                feynkit_graph::symbols::ubar(),
                feynkit_graph::symbols::v(),
                feynkit_graph::symbols::vbar(),
            ] {
                assert_eq!(function_count(diagram.projector(), wavefunction), 1);
            }
            assert_eq!(
                function_count(diagram.projector(), feynkit_graph::symbols::hedge_index()),
                4
            );
        }

        let explicit = generator
            .generate(
                &process,
                &options.projector(Atom::one()),
                GenerationType::Amplitude,
            )
            .unwrap();
        assert!(
            explicit
                .diagrams
                .iter()
                .all(|diagram| diagram.projector() == &Atom::one())
        );
    }

    #[test]
    fn compton_cuts_preserve_the_requested_particle_charge() {
        let generator = Generator::new(fermion_model());
        for pdg in [1_i64, -1] {
            let process = Process::new([pdg, 22], [pdg, 22]);
            let generated = generator
                .generate(
                    &process,
                    &GenerationOptions::default()
                        .with_loop_count(1, 1)
                        .unwrap()
                        .threads(1)
                        .max_vertices(4),
                    GenerationType::CrossSection,
                )
                .unwrap();
            assert!(!generated.diagrams.is_empty());
            for diagram in generated.diagrams {
                assert!(!diagram.cuts().is_empty());
                for cut in diagram.cuts() {
                    let mut particles = diagram
                        .cut_particles(cut)
                        .unwrap()
                        .into_iter()
                        .map(|id| diagram.model().particle_by_id(id).unwrap().pdg_code)
                        .collect::<Vec<_>>();
                    particles.sort();
                    assert_eq!(particles, vec![pdg, 22]);
                    // Reversing a stored edge conjugates its species and its
                    // cut endpoint together, preserving the physical state.
                    let reversed = diagram.reverse_edge(cut.cut[0].edge).unwrap();
                    let mut reversed_particles = reversed
                        .cut_particles(&reversed.cuts()[0])
                        .unwrap()
                        .into_iter()
                        .map(|id| reversed.model().particle_by_id(id).unwrap().pdg_code)
                        .collect::<Vec<_>>();
                    reversed_particles.sort();
                    assert_eq!(reversed_particles, particles);
                }
            }
        }
    }

    #[test]
    fn sewing_a_single_external_fermion_chain_has_even_ordering() {
        let model = fermion_model();
        let generator = Generator::new(model);
        for pdg in [1_i64, -1] {
            let process = Process::new([pdg, 22], [pdg, 22]);
            let generated = generator
                .generate(
                    &process,
                    &GenerationOptions::default()
                        .with_loop_count(1, 1)
                        .unwrap()
                        .threads(1)
                        .max_vertices(4),
                    GenerationType::CrossSection,
                )
                .unwrap();
            assert!(!generated.diagrams.is_empty());
            for diagram in generated.diagrams {
                // One external chain closes into one cycle without exchanging
                // any external fermions. Its ordering sign must be positive.
                let positive = symbolica::function!(
                    symbol!("feynkit_generator_factor::ExternalFermionOrderingSign"),
                    Atom::num(1)
                );
                assert!(diagram.overall_factor().contains(positive.as_view()));
            }
        }
    }

    #[test]
    fn normalized_cross_section_fermion_chains_keep_valid_cuts() {
        let generator = Generator::new(fermion_model());
        let process = Process::new([1_i64, -1], [22_i64]);
        let generated = generator
            .generate(
                &process,
                &GenerationOptions::default()
                    .with_loop_count(0, 1)
                    .unwrap()
                    .threads(1)
                    .max_vertices(4),
                GenerationType::CrossSection,
            )
            .unwrap();
        assert!(!generated.diagrams.is_empty());

        for diagram in &generated.diagrams {
            diagram.validate().unwrap();
            assert!(!diagram.cuts().is_empty());

            let external_edges = diagram
                .edges()
                .filter_map(|(edge, _, data)| data.external.is_some().then_some(edge))
                .collect::<BTreeSet<_>>();
            for candidate in diagram.topology_threshold_candidates() {
                assert!(
                    candidate
                        .cut
                        .iter()
                        .all(|half_edge| !external_edges.contains(&half_edge.edge)),
                    "topology thresholds must be cuts of the sewn internal attachment nodes"
                );
            }

            let fermion_edges = diagram
                .edges()
                .filter(|(_, _, edge)| edge.external.is_none())
                .filter_map(|(edge_id, endpoints, edge)| {
                    diagram
                        .model()
                        .particle_by_id(edge.particle)
                        .unwrap()
                        .is_fermion()
                        .then_some((
                            edge_id,
                            endpoints.source.expect("sewn source").0,
                            endpoints.target.expect("sewn target").0,
                            edge.particle,
                        ))
                })
                .collect::<Vec<_>>();
            let mut remaining = fermion_edges
                .iter()
                .map(|(edge, _, _, _)| *edge)
                .collect::<BTreeSet<_>>();
            // Sewing preserves the ordinary particle flow at both attachment
            // vertices, including initial-state carriers.
            for (node, _, _) in diagram.underlying().iter_nodes() {
                let signed_pdgs = diagram
                    .underlying()
                    .iter_crown(node)
                    .map(|hedge| {
                        let edge = &diagram.underlying()[diagram.underlying()[&hedge]];
                        let particle = diagram.model().particle_by_id(edge.particle).unwrap();
                        let pdg = particle.pdg_code;
                        if !particle.is_fermion() {
                            0
                        } else if diagram.underlying().flow(hedge) == Flow::Source {
                            pdg
                        } else {
                            -pdg
                        }
                    })
                    .sum::<i64>();
                assert_eq!(
                    signed_pdgs, 0,
                    "fermion flow must be conserved at each interaction"
                );
            }
            while let Some(start) = remaining.iter().next().copied() {
                let (_, source, target, _) = fermion_edges
                    .iter()
                    .find(|(edge, _, _, _)| *edge == start)
                    .copied()
                    .unwrap();
                let mut stack = vec![source, target];
                let mut component_particles = BTreeSet::new();
                while let Some(vertex) = stack.pop() {
                    for (edge, source, target, particle) in &fermion_edges {
                        if remaining.contains(edge) && (*source == vertex || *target == vertex) {
                            remaining.remove(edge);
                            component_particles.insert(*particle);
                            stack.push(*source);
                            stack.push(*target);
                        }
                    }
                }
                assert_eq!(
                    component_particles.len(),
                    1,
                    "every finalized fermion chain must use one particle species"
                );
            }
        }
    }

    #[test]
    fn cut_partitions_use_sewn_attachment_nodes() {
        let model = fermion_model();
        let resolved = ResolvedProcess::new(
            &model,
            &Process::new([22_i64, 22], [22_i64, 22]),
            GenerationType::CrossSection,
            &GenerationOptions::default(),
        )
        .unwrap();
        let mut graph = Graph::new();
        let external = |index, state| {
            ColoredNode::External(ExternalNode {
                index,
                state,
                particle: particle(&model, 22),
            })
        };
        let incoming_0 = graph.add_node(external(0, ExternalState::Incoming));
        let incoming_1 = graph.add_node(external(1, ExternalState::Incoming));
        let outgoing_0 = graph.add_node(external(2, ExternalState::Outgoing));
        let outgoing_1 = graph.add_node(external(3, ExternalState::Outgoing));
        let left_0 = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let left_1 = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let right_0 = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let right_1 = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let photon = edge(&model, 22);
        let external_edges = BTreeSet::from([
            EdgeId(graph.add_edge(incoming_0, left_0, false, photon).unwrap()),
            EdgeId(graph.add_edge(incoming_1, left_1, false, photon).unwrap()),
            EdgeId(graph.add_edge(right_0, outgoing_0, false, photon).unwrap()),
            EdgeId(graph.add_edge(right_1, outgoing_1, false, photon).unwrap()),
        ]);
        let bridge_0 = EdgeId(graph.add_edge(left_0, right_0, false, photon).unwrap());
        let bridge_1 = EdgeId(graph.add_edge(left_1, right_1, false, photon).unwrap());
        graph.add_edge(left_0, left_1, false, photon).unwrap();
        graph.add_edge(right_0, right_1, false, photon).unwrap();

        // The former pre-sewing convention fixed only the degree-one external
        // vertices. It therefore admitted partitions which crossed an initial
        // external edge and two or more internal edges; stripping that edge in
        // the runtime bridge turned them into non-s-t threshold candidates.
        let old_graph = graph.build().with_external_vertices();
        let mut old_source = Vec::new();
        let mut old_target = Vec::new();
        for (node, _, color) in old_graph.iter_nodes() {
            let PlainVertex::External(external) = color else {
                continue;
            };
            match external.state {
                ExternalState::Incoming => old_source.push(node),
                ExternalState::Outgoing => old_target.push(node),
            }
        }
        assert!(
            old_graph
                .all_cuts_from_ids(&old_source, &old_target)
                .into_iter()
                .any(|(_, cut, _)| {
                    let crossing = cut
                        .left
                        .included_iter()
                        .map(|hedge| EdgeId(old_graph[&hedge].0))
                        .collect::<BTreeSet<_>>();
                    crossing.iter().any(|edge| external_edges.contains(edge))
                        && crossing
                            .iter()
                            .filter(|edge| !external_edges.contains(edge))
                            .count()
                            > 1
                })
        );

        let inventory = resolved
            .cut_partitions(&graph.build(), &model, &GenerationOptions::default())
            .unwrap();
        assert_eq!(inventory.topology_threshold_candidates.len(), 1);
        let crossing = inventory.topology_threshold_candidates[0]
            .cut
            .iter()
            .map(|half_edge| half_edge.edge)
            .collect::<BTreeSet<_>>();
        assert_eq!(crossing, BTreeSet::from([bridge_0, bridge_1]));
        assert!(crossing.is_disjoint(&external_edges));
    }

    #[test]
    fn default_loop_routing_uses_canonical_topology_order() {
        let model = Arc::new(scalar_model());
        let particle = model.particle_id("phi").unwrap();
        let interaction = model.vertex_rule_id("V1").unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "routing-order");
        let vertices = (0..4)
            .map(|vertex| {
                builder.add_vertex(DiagramVertex::interaction(
                    format!("v{vertex}"),
                    interaction,
                ))
            })
            .collect::<Vec<_>>();
        for (edge, source, target) in [(0, 2, 3), (1, 1, 3), (2, 1, 2), (3, 0, 3), (4, 0, 1)] {
            assert_eq!(
                builder
                    .add_edge(
                        vertices[source],
                        vertices[target],
                        DiagramEdge::new(particle, false),
                    )
                    .unwrap(),
                EdgeId(edge),
            );
        }
        let stable_id_default = builder.build().unwrap();
        assert_eq!(
            stable_id_default.loop_momentum_basis().loop_edges,
            vec![EdgeId(0), EdgeId(1)],
        );

        let canonical = stable_id_default
            .with_canonical_loop_momentum_basis()
            .unwrap();
        assert_eq!(
            canonical.loop_momentum_basis().tree_edges,
            vec![EdgeId(0), EdgeId(2), EdgeId(3)],
        );
        assert_eq!(
            canonical.loop_momentum_basis().loop_edges,
            vec![EdgeId(1), EdgeId(4)],
        );
    }

    #[test]
    fn default_loop_routing_matches_gl_multichannel_traversal() {
        let model = Arc::new(scalar_model());
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "gl-like-routing");
        let vertices = (0..8)
            .map(|vertex| {
                builder.add_vertex(DiagramVertex {
                    name: format!("v{vertex}"),
                    interaction: None,

                    numerator: Atom::one(),
                })
            })
            .collect::<Vec<_>>();
        for index in 0..2 {
            assert_eq!(
                builder
                    .add_generation_external(
                        vertices[7],
                        Flow::Sink,
                        DiagramEdge::new(particle, false),
                        VertexSlot(0),
                        ExternalLeg {
                            name: format!("out{index}"),
                            index,
                            state: ExternalState::Outgoing,
                            connection: index,
                        },
                    )
                    .unwrap(),
                EdgeId(index),
            );
        }
        for (edge, source, target) in [
            (0, 0, 1),
            (1, 0, 3),
            (2, 5, 0),
            (3, 2, 1),
            (4, 1, 5),
            (5, 4, 2),
            (6, 2, 7),
            (7, 3, 4),
            (8, 3, 6),
            (9, 4, 5),
        ] {
            assert_eq!(
                builder
                    .add_edge(
                        vertices[source],
                        vertices[target],
                        DiagramEdge::new(particle, false),
                    )
                    .unwrap(),
                EdgeId(edge + 2),
            );
        }

        let routed = builder.canonical_loop_momentum_basis().build().unwrap();

        assert_eq!(
            routed.loop_momentum_basis().loop_edges,
            vec![EdgeId(4), EdgeId(5), EdgeId(9)],
        );
    }

    /// Cross-section integrands rescale the loop momenta about the origin of the
    /// default basis, which lies inside every cut energy surface only when no
    /// cut edge carries a partial beam momentum such as `K + P1`.
    #[test]
    fn default_cross_section_routing_keeps_partial_beam_momenta_off_cut_edges() {
        let model = standard_model();
        let qcd = ["g", "t", "t~", "ghG", "ghG~"];
        let qed = ["e-", "e+", "a"];
        for (incoming, outgoing, selection, loops) in [
            (["a", "a"], ["t", "t~"], &["a", "t", "t~"][..], 1),
            (["Z", "Z"], ["t", "t~"], &["Z", "t", "t~", "H"][..], 1),
            (["g", "g"], ["t", "t~"], &qcd[..], 1),
            (["e-", "e-"], ["e-", "e-"], &qed[..], 1),
            // Two loops add real-emission diagrams with several cuts.
            (["e-", "e-"], ["e-", "e-"], &qed[..], 2),
            (["e+", "e-"], ["a", "a"], &qed[..], 2),
            (["a", "a"], ["t", "t~"], &["a", "t", "t~", "g"][..], 2),
        ] {
            let process = Process::new(incoming, outgoing)
                .with_particle_filters(
                    None,
                    Some(selection.iter().map(|&name| name.into()).collect()),
                    None,
                    Vec::new(),
                )
                .unwrap();
            for symmetrize in [false, true] {
                let options = GenerationOptions::default()
                    .with_loop_count(loops, loops)
                    .unwrap()
                    .threads(1)
                    .symmetrize_initial(symmetrize);
                let result = process
                    .generate_cross_section(model.clone(), &options)
                    .unwrap();
                assert!(!result.diagrams.is_empty(), "{incoming:?} > {outgoing:?}");
                let several_cuts = result
                    .diagrams
                    .iter()
                    .any(|diagram| diagram.cuts().len() > 1);
                assert_eq!(several_cuts, loops > 1, "{incoming:?} > {outgoing:?}");
                for diagram in &result.diagrams {
                    for cut in diagram.cuts() {
                        for half_edge in &cut.cut {
                            let momentum =
                                &diagram.loop_momentum_basis().edge_signatures[&half_edge.edge];
                            let shift = &momentum.external;
                            assert!(
                                shift.iter().all(|sign| Some(sign) == shift.get(0)),
                                "{incoming:?} > {outgoing:?} L{loops} {:?}: cut edge {} carries {}",
                                diagram.name(),
                                half_edge.edge.0,
                                momentum.format_momentum(),
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn finalizes_vector_projectors_and_leaves_scalar_states_trivial() {
        let options = GenerationOptions::default()
            .with_loop_count(1, 1)
            .unwrap()
            .threads(1);
        let vectors = Generator::new(fermion_model())
            .generate(&photon_process(), &options, GenerationType::Amplitude)
            .unwrap();
        assert!(!vectors.diagrams.is_empty());
        for diagram in &vectors.diagrams {
            for (head, count) in [
                (feynkit_graph::symbols::epsilon(), 1),
                (feynkit_graph::symbols::epsilonbar(), 1),
                (feynkit_graph::symbols::hedge_index(), 2),
            ] {
                assert_eq!(function_count(diagram.projector(), head), count);
            }
        }

        let scalars = Generator::new(scalar_model())
            .generate(
                &Process::new(["phi"], ["phi"]),
                &options,
                GenerationType::Amplitude,
            )
            .unwrap();
        assert!(!scalars.diagrams.is_empty());
        assert!(
            scalars
                .diagrams
                .iter()
                .all(|diagram| diagram.projector() == &Atom::one())
        );
    }

    #[test]
    fn forced_cuts_retain_only_the_requested_typed_edge_set() {
        let generator = Generator::new(scalar_model());
        let process = scalar_cross_section_process();
        let options = GenerationOptions::default()
            .with_loop_count(1, 1)
            .unwrap()
            .threads(1)
            .max_vertices(4);
        let generated = generator
            .generate(&process, &options, GenerationType::CrossSection)
            .unwrap();
        let requested: BTreeSet<EdgeId> = generated
            .diagrams
            .iter()
            .flat_map(|diagram| diagram.cuts())
            .next()
            .unwrap()
            .cut
            .iter()
            .map(|half_edge| half_edge.edge)
            .collect();

        let forced = generator
            .generate(
                &process,
                &options.clone().forced_cuts([requested.clone()]),
                GenerationType::CrossSection,
            )
            .unwrap();

        assert!(!forced.diagrams.is_empty());
        assert!(forced.diagrams.iter().all(|diagram| {
            !diagram.cuts().is_empty()
                && diagram.cuts().iter().all(|cut| {
                    cut.cut
                        .iter()
                        .map(|half_edge| half_edge.edge)
                        .collect::<BTreeSet<_>>()
                        == requested
                })
        }));
    }

    #[test]
    fn unknown_forced_cut_reports_the_typed_edge_set() {
        let unknown = BTreeSet::from([EdgeId(usize::MAX)]);
        let result = Generator::new(scalar_model()).generate(
            &scalar_cross_section_process(),
            &GenerationOptions::default()
                .with_loop_count(1, 1)
                .unwrap()
                .threads(1)
                .max_vertices(4)
                .forced_cuts([unknown]),
            GenerationType::CrossSection,
        );

        assert!(matches!(
            result,
            Err(GenerationError::UnknownForcedCut { edges })
                if edges == vec![usize::MAX]
        ));
    }

    #[test]
    fn forced_cut_filtering_precedes_numerator_grouping() {
        let generator = Generator::new(scalar_model());
        let process = scalar_cross_section_process();
        let options = GenerationOptions::default()
            .with_loop_count(1, 1)
            .unwrap()
            .threads(1)
            .max_vertices(4);
        let generated = generator
            .generate(&process, &options, GenerationType::CrossSection)
            .unwrap();
        let cut_sets = generated
            .diagrams
            .iter()
            .flat_map(|diagram| diagram.cuts())
            .map(|cut| {
                cut.cut
                    .iter()
                    .map(|half_edge| half_edge.edge)
                    .collect::<BTreeSet<_>>()
            })
            .collect::<BTreeSet<_>>();
        assert!(
            cut_sets.len() > 1,
            "fixture must exercise actual cut filtering"
        );
        let requested = cut_sets.into_iter().next().unwrap();

        let grouped = generator
            .generate(
                &process,
                &options.forced_cuts([requested.clone()]).numerator_grouping(
                    NumeratorGrouping::Identical(GraphGroupingOptions::default()),
                ),
                GenerationType::CrossSection,
            )
            .unwrap();

        grouped.validate_groups().unwrap();
        assert!(!grouped.groups.is_empty());
        assert!(grouped.diagrams.iter().all(|diagram| {
            !diagram.cuts().is_empty()
                && diagram.cuts().iter().all(|cut| {
                    cut.cut
                        .iter()
                        .map(|half_edge| half_edge.edge)
                        .collect::<BTreeSet<_>>()
                        == requested
                })
        }));
    }

    #[test]
    fn rejects_model_bound_process_and_filter_selectors_from_another_model() {
        let source = scalar_model();
        let mut target_json: serde_json::Value =
            serde_json::from_str(&source.to_json().unwrap()).unwrap();
        target_json["name"] = serde_json::Value::String("other_phi3".to_owned());
        let target = Model::from_json(&serde_json::to_string(&target_json).unwrap()).unwrap();
        assert_ne!(source.fingerprint(), target.fingerprint());

        let foreign_particle =
            ParticleSelector::by_id(&source, source.particle_id("phi").unwrap()).unwrap();
        let process = Process::new([foreign_particle.clone()], ["phi"]);
        assert!(matches!(
            Generator::new(target.clone()).generate(
                &process,
                &GenerationOptions::default(),
                GenerationType::Amplitude
            ),
            Err(GenerationError::Selector(SelectorError::ModelMismatch {
                kind: "particle",
                ..
            }))
        ));

        let foreign_vertex =
            VertexSelector::by_id(&source, source.vertex_rule_id("V1").unwrap()).unwrap();
        for filter in [
            GenerationFilter::ParticleVeto(vec![foreign_particle]),
            GenerationFilter::VertexAllow(vec![foreign_vertex]),
        ] {
            let options = GenerationOptions::default().with_graph_filter(filter);
            assert!(matches!(
                Generator::new(target.clone()).generate(
                    &Process::new(["phi"], ["phi", "phi"]),
                    &options,
                    GenerationType::Amplitude
                ),
                Err(GenerationError::Selector(
                    SelectorError::ModelMismatch { .. }
                ))
            ));
        }
    }

    #[test]
    fn empty_generation_returns_a_completed_result() {
        let generator = Generator::new(scalar_model());
        let process = Process::new(["phi"], ["phi", "phi"]);
        for options in [
            GenerationOptions::default().max_vertices(0),
            GenerationOptions::default().with_graph_filter(GenerationFilter::ParticleVeto(vec![
                ParticleSelector::Name("phi".to_owned()),
            ])),
        ] {
            let result = generator
                .generate(&process, &options, GenerationType::Amplitude)
                .unwrap();
            assert!(result.report.completed);
            assert!(result.diagrams.is_empty());
            assert!(result.groups.is_empty());
            assert_eq!(result.report.topology_count, 0);
            assert_eq!(result.report.interaction_assignment_count, 0);
            assert_eq!(result.report.retained_count, 0);
        }
    }

    #[test]
    fn cancelled_generation_returns_an_incomplete_result() {
        let cancellation = CancellationToken::new();
        cancellation.cancel();
        let generated = Generator::new(scalar_model())
            .generate(
                &Process::new(["phi"], ["phi", "phi"]),
                &GenerationOptions::default().cancellation_token(cancellation),
                GenerationType::Amplitude,
            )
            .unwrap();

        assert!(!generated.report.completed);
        assert!(generated.diagrams.is_empty());
    }

    /// The overall factor with the generator's symbolic factors replaced by their arguments.
    fn numeric_factor(diagram: &FeynmanDiagram) -> Atom {
        diagram
            .overall_factor()
            .replace_map(|term, _, out| {
                if let AtomView::Fun(function) = term
                    && function.get_symbol().get_namespace() == "feynkit_generator_factor"
                {
                    **out = function.iter().next().unwrap().to_owned();
                }
            })
            .expand()
    }

    #[test]
    fn split_vertices_with_identical_fields_give_distinct_diagrams() {
        let mut definition = serde_json::to_value(scalar_model()).unwrap();
        let rules = definition["vertex_rules"].as_array_mut().unwrap();
        let mut counterterm = rules[0].clone();
        counterterm["name"] = serde_json::json!("V1_CT");
        rules.push(counterterm);
        let model: Model = serde_json::from_value(definition).unwrap();
        let generated = Generator::new(model.clone())
            .generate(
                &Process::new(["phi"], ["phi"]),
                &GenerationOptions::default()
                    .with_loop_count(1, 1)
                    .unwrap()
                    .threads(1),
                GenerationType::Amplitude,
            )
            .unwrap();
        // The bubble with either rule at either vertex: the legs distinguish the vertices, and
        // the two internal lines are interchangeable.
        let mut rules = generated
            .diagrams
            .iter()
            .map(|diagram| {
                assert_eq!(diagram.symmetry_factor(), 2);
                assert_eq!(numeric_factor(diagram), Atom::num(1) / Atom::num(2));
                diagram
                    .vertices()
                    .map(|(_, vertex)| {
                        model
                            .vertex_rule_by_id(vertex.interaction.unwrap())
                            .unwrap()
                            .name
                            .clone()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        for rules in &mut rules {
            rules.sort();
        }
        rules.sort();
        assert_eq!(
            rules,
            [
                ["V1", "V1"],
                ["V1", "V1_CT"],
                ["V1", "V1_CT"],
                ["V1_CT", "V1_CT"]
            ]
        );
    }

    #[test]
    fn symmetrized_external_momenta_are_assigned_canonically() {
        let generator = Generator::new(Model::qcd());
        let process = Process::new(["g", "g"], ["g", "g"]);
        let options = GenerationOptions::default()
            .with_loop_count(0, 0)
            .unwrap()
            .symmetrize_initial(true)
            .symmetrize_final(true);
        let ids = |threads| {
            generator
                .generate(
                    &process,
                    &options.clone().threads(threads),
                    GenerationType::Amplitude,
                )
                .unwrap()
                .diagrams
                .iter()
                .map(|diagram| {
                    let mut legs = diagram
                        .edges()
                        .filter_map(|(_, _, edge)| edge.external.as_ref().map(|leg| leg.index))
                        .collect::<Vec<_>>();
                    legs.sort();
                    assert_eq!(legs, [0, 1, 2, 3]);
                    diagram.id().unwrap()
                })
                .collect::<Vec<_>>()
        };
        let serial = ids(1);
        // s channel, t channel (the u channel exchanges the final gluons) and the contact vertex
        assert_eq!(serial.len(), 3);
        assert_eq!(serial, ids(4));
    }

    #[test]
    fn output_metadata_and_explicit_selection_are_independent() {
        let model = Model::qed();
        let process = Process::new(["e-", "e+"], ["a", "a"]);
        let minimal = GenerationOptions::default().output(GenerationOutput::Minimal);
        let all = process.generate_diagrams(model.clone(), &minimal).unwrap();
        assert_eq!(all.diagrams.len(), 2);
        assert!(all.diagrams.iter().all(|diagram| diagram.name().is_none()));
        assert!(
            all.groups
                .iter()
                .flat_map(|group| &group.members)
                .all(|member| member.source_identity.is_none())
        );
        let chosen = all.diagrams[0].id().unwrap();
        let selected = process
            .generate_diagrams(model.clone(), &minimal.clone().select_diagram_ids([chosen]))
            .unwrap();
        assert_eq!(selected.diagrams.len(), 1);
        assert_eq!(selected.diagrams[0].id().unwrap(), chosen);
        assert!(selected.diagrams[0].name().is_none());
        let vetoed = process
            .generate_diagrams(model.clone(), &minimal.clone().veto_diagram_ids([chosen]))
            .unwrap();
        assert_eq!(vetoed.diagrams.len(), 1);
        assert_ne!(vetoed.diagrams[0].id().unwrap(), chosen);
        for empty in [
            minimal.clone().select_diagram_ids([]),
            minimal.clone().select_diagrams(Vec::<String>::new()),
        ] {
            assert!(
                process
                    .generate_diagrams(model.clone(), &empty)
                    .unwrap()
                    .diagrams
                    .is_empty()
            );
        }
        let named = minimal.clone().output(GenerationOutput::Named {
            prefix: "test_".into(),
        });
        let selected = process
            .generate_diagrams(model.clone(), &named.clone().select_diagrams(["test_0"]))
            .unwrap();
        assert_eq!(selected.diagrams.len(), 1);
        assert_eq!(selected.diagrams[0].name(), Some("test_0"));
        assert_eq!(selected.diagrams[0].id().unwrap(), chosen);
        // Both requests match: ID routing takes precedence over the invalid name map.
        let routed = named
            .with_loop_momentum_basis(chosen, [])
            .with_named_loop_momentum_basis("test_0", [EdgeId(usize::MAX)]);
        assert_eq!(
            process
                .generate_diagrams(model, &routed)
                .unwrap()
                .diagrams
                .len(),
            2
        );
    }

    #[test]
    fn interned_numerator_colors_preserve_old_string_order_and_collisions() {
        // No interactions are needed to exercise the color alphabet. Include a
        // particle name equal to another role's mass/spin spelling and spins
        // whose decimal order differs from their numeric order.
        let mut definition: serde_json::Value =
            serde_json::from_str(&scalar_model().to_json().unwrap()).unwrap();
        for field in [
            "vertex_rules",
            "propagators",
            "lorentz_structures",
            "couplings",
        ] {
            definition[field] = serde_json::json!([]);
        }
        let template = definition["particles"][0].clone();
        definition["particles"] = serde_json::json!([]);
        for (index, (name, spin)) in [("M | 2", 1), ("ten", 10), ("two", 2)]
            .into_iter()
            .enumerate()
        {
            let mut particle = template.clone();
            particle["name"] = name.into();
            particle["antiname"] = name.into();
            particle["spin"] = spin.into();
            particle["pdg_code"] = (100 + index).into();
            definition["particles"]
                .as_array_mut()
                .unwrap()
                .push(particle);
        }
        let collision_model = Model::from_json(&definition.to_string()).unwrap();
        for model in [collision_model, Model::qed(), Model::qcd()] {
            let colors = NumeratorColors::new(&model).unwrap();
            let mut observations = Vec::new();
            for (index, particle) in model.particles().iter().enumerate() {
                let id = model.particle_id_at(index).unwrap();
                let external = Endpoint::External(ExternalNode {
                    index: 0,
                    state: ExternalState::Outgoing,
                    particle: id,
                });
                for (source, sink) in [
                    (Endpoint::Vertex(0), Endpoint::Vertex(1)),
                    (external, Endpoint::Vertex(0)),
                    (Endpoint::Vertex(0), external),
                ] {
                    let edge = GeneratedEdge {
                        source,
                        sink,
                        particle: id,
                        directed: !model.particle_is_self_conjugate(id),
                    };
                    for exact in [false, true] {
                        let entering = if edge.is_external() && edge.sink.external().is_some() {
                            model.particle_by_id(particle.antiparticle).unwrap()
                        } else {
                            particle
                        };
                        let old = if edge.is_external() || exact {
                            entering.name.clone()
                        } else {
                            format!(
                                "{} | {}",
                                model.parameter_by_id(particle.mass).unwrap().name,
                                particle.spin
                            )
                        };
                        observations.push((colors.edge(&edge, &model, exact).unwrap(), old));
                    }
                }
            }
            for (left_rank, left_text) in &observations {
                for (right_rank, right_text) in &observations {
                    assert_eq!(
                        left_rank.cmp(right_rank),
                        left_text.cmp(right_text),
                        "{left_text:?} vs {right_text:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn thread_counts_preserve_diagrams_and_progress_order() {
        let generator = Generator::new(Model::qcd());
        let process = Process::new(["g"], ["g"])
            .with_particle_filters(
                Some(
                    ["u", "d", "c", "s", "t"]
                        .into_iter()
                        .map(Into::into)
                        .collect(),
                ),
                None,
                None,
                vec![],
            )
            .unwrap();
        let options = GenerationOptions::default()
            .with_loop_count(2, 2)
            .unwrap()
            .with_graph_filter(GenerationFilter::MaxNumberOfBridges(0));
        let serial = generator
            .generate(
                &process,
                &options.clone().threads(1),
                GenerationType::Amplitude,
            )
            .unwrap();
        let progress = Arc::new(Mutex::new(Vec::new()));
        let progress_sink = Arc::clone(&progress);
        let parallel = generator
            .generate(
                &process,
                &options.threads(4).progress(move |update| {
                    if update.stage == "numerators" {
                        progress_sink.lock().unwrap().push(update.completed);
                    }
                    crate::GenerationControl::Continue
                }),
                GenerationType::Amplitude,
            )
            .unwrap();
        assert!(!serial.diagrams.is_empty());
        assert_eq!(serial.report, parallel.report);
        assert_eq!(
            *progress.lock().unwrap(),
            (0..=parallel.diagrams.len()).collect::<Vec<_>>()
        );
        for (serial, parallel) in serial.diagrams.iter().zip(&parallel.diagrams) {
            assert_eq!(serial.to_json().unwrap(), parallel.to_json().unwrap());
        }
        assert_eq!(serial.groups.len(), parallel.groups.len());
        for (serial, parallel) in serial.groups.iter().zip(&parallel.groups) {
            assert_eq!(serial.master, parallel.master);
            assert_eq!(serial.members.len(), parallel.members.len());
            for (serial, parallel) in serial.members.iter().zip(&parallel.members) {
                assert_eq!(serial.source_identity, parallel.source_identity);
                assert_eq!(serial.source_diagram, parallel.source_diagram);
                assert_eq!(serial.ratio, parallel.ratio);
                assert_eq!(serial.overall_factor, parallel.overall_factor);
            }
        }
    }

    #[test]
    fn zero_color_filter_is_opt_in_for_generated_qcd_tadpoles() {
        let model = Arc::new(Model::qcd());
        let process = Process::new([] as [&str; 0], ["g"]);
        let options = GenerationOptions::default()
            .with_loop_count(1, 1)
            .unwrap()
            .max_vertices(1)
            .allow_self_loops(true)
            .allow_zero_flow_edges(true)
            .threads(1);
        let unfiltered = process.generate_diagrams(model.clone(), &options).unwrap();
        assert!(unfiltered.report.completed);
        assert!(!unfiltered.diagrams.is_empty());
        assert_eq!(unfiltered.report.zero_numerator_count, 0);
        assert!(
            unfiltered
                .diagrams
                .iter()
                .any(|diagram| !diagram.numerator().is_zero())
        );

        let updates = Arc::new(Mutex::new(Vec::new()));
        let captured = updates.clone();
        let filtered_options = options
            .clone()
            .filter_zero_color(true)
            .progress(move |update| {
                if update.stage == "filter_zero_color" {
                    captured.lock().unwrap().push(update);
                }
                crate::GenerationControl::Continue
            });
        let filtered = process
            .generate_diagrams(model.clone(), &filtered_options)
            .unwrap();
        assert!(filtered.report.completed);
        assert!(filtered.diagrams.is_empty());
        assert_eq!(
            filtered.report.zero_numerator_count,
            unfiltered.diagrams.len()
        );
        filtered.validate_groups().unwrap();
        let updates = updates.lock().unwrap();
        assert_eq!(
            updates
                .iter()
                .map(|update| (update.completed, update.zero_numerator_count))
                .collect::<Vec<_>>(),
            (0..=unfiltered.diagrams.len())
                .map(|count| (count, count))
                .collect::<Vec<_>>()
        );
        assert!(
            updates
                .iter()
                .all(|update| update.total == Some(unfiltered.diagrams.len()))
        );

        let amplitude = process
            .generate_amplitude(model.clone(), &options, Default::default())
            .unwrap();
        assert_eq!(amplitude.diagrams().len(), unfiltered.diagrams.len());
        assert!(matches!(
            process.generate_amplitude(model, &options.filter_zero_color(true), Default::default()),
            Err(GenerationError::Amplitude(error)) if matches!(*error, feynkit_amplitude::AmplitudeError::Empty)
        ));
    }

    #[test]
    fn zero_color_filter_progress_can_cancel_algebra() {
        let process = Process::new([] as [&str; 0], ["g"]);
        let result = process
            .generate_diagrams(
                Model::qcd(),
                &GenerationOptions::default()
                    .with_loop_count(1, 1)
                    .unwrap()
                    .max_vertices(1)
                    .allow_self_loops(true)
                    .allow_zero_flow_edges(true)
                    .threads(1)
                    .filter_zero_color(true)
                    .progress(|update| {
                        if update.stage == "filter_zero_color" && update.completed == 1 {
                            crate::GenerationControl::Cancel
                        } else {
                            crate::GenerationControl::Continue
                        }
                    }),
            )
            .unwrap();
        assert!(!result.report.completed);
        assert_eq!(result.report.zero_numerator_count, 1);
        assert!(!result.diagrams.is_empty());
        result.validate_groups().unwrap();
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn color_filter_uses_the_requested_generation_pool() {
        for threads in [1, 4] {
            let observed = Arc::new(Mutex::new(Vec::new()));
            let sink = observed.clone();
            let options = GenerationOptions::default()
                .with_loop_count(1, 1)
                .unwrap()
                .threads(threads)
                .filter_zero_color(true)
                .progress(move |progress| {
                    if progress.stage == "filter_zero_color" && progress.completed > 0 {
                        sink.lock().unwrap().push(rayon::current_num_threads());
                    }
                    crate::GenerationControl::Continue
                });
            let generated = Generator::new(Model::yang_mills())
                .generate(
                    &Process::new(["g"], ["g"]),
                    &options,
                    GenerationType::Amplitude,
                )
                .unwrap();
            assert!(generated.report.completed);
            let observed = observed.lock().unwrap();
            assert!(!observed.is_empty());
            assert!(observed.iter().all(|&actual| actual == threads));
        }
    }

    #[test]
    fn numerator_progress_can_cancel_construction() {
        let generated = Generator::new(Model::yang_mills())
            .generate(
                &Process::new(["g"], ["g"]),
                &GenerationOptions::default()
                    .with_loop_count(2, 2)
                    .unwrap()
                    .threads(4)
                    .progress(|update| {
                        if update.stage == "numerators" && update.completed == 2 {
                            crate::GenerationControl::Cancel
                        } else {
                            crate::GenerationControl::Continue
                        }
                    }),
                GenerationType::Amplitude,
            )
            .unwrap();
        assert!(!generated.report.completed);
        assert!(!generated.diagrams.is_empty());
        assert!(generated.report.retained_count < generated.report.interaction_assignment_count);
        generated.validate_groups().unwrap();
    }

    #[test]
    fn cancellation_callback_returns_an_incomplete_result() {
        let generated = Generator::new(scalar_model())
            .generate(
                &Process::new(["phi"], ["phi", "phi"]),
                &GenerationOptions::default().cancellation_check(|| true),
                GenerationType::Amplitude,
            )
            .unwrap();

        assert!(!generated.report.completed);
        assert!(generated.diagrams.is_empty());
    }

    #[test]
    fn vertex_signatures_present_ufo_particles_as_incoming_stubs() {
        let model = fermion_model();
        let options =
            GenerationOptions::default().with_graph_filter(GenerationFilter::VertexAllow(vec![
                VertexSelector::Name("V".to_owned()),
            ]));
        let options = options.resolve_selectors(&model).unwrap();
        let signatures = VertexSignatures::new(&model, &options).unwrap();
        assert_eq!(signatures.rules, [vertex(&model, "V")]);
        // UFO particles enter the vertex: the antifermion f~ is an outgoing f.
        assert_eq!(
            signatures.stubs,
            [vec![
                ParticleStub::outgoing(particle(&model, 1)),
                ParticleStub::incoming(particle(&model, 1)),
                ParticleStub::undirected(particle(&model, 22)),
            ]]
        );
    }

    #[test]
    fn external_fermion_chains_are_oriented_consistently() {
        let model = fermion_model();
        let process = ResolvedProcess::new(
            &model,
            &Process::new([1_i64, -1], Vec::<i64>::new()),
            GenerationType::Amplitude,
            &GenerationOptions::default(),
        )
        .unwrap();
        let mut graph = Graph::new();
        let particle_vertex = graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Incoming,
            particle: particle(&model, 1),
        }));
        let left = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let right = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let antiparticle = graph.add_node(ColoredNode::External(ExternalNode {
            index: 1,
            state: ExternalState::Incoming,
            particle: particle(&model, -1),
        }));
        let fermion = particle(&model, 1);
        graph
            .add_edge(particle_vertex, left, true, fermion)
            .unwrap();
        graph
            .add_edge(right, left, true, particle(&model, -1))
            .unwrap();
        graph.add_edge(right, antiparticle, true, fermion).unwrap();

        let normalized = process
            .normalize_fermion_flows(&graph.build(), &model)
            .unwrap();
        assert_eq!(graph.ends(&normalized.graph, 0), (particle_vertex, left));
        assert_eq!(graph.ends(&normalized.graph, 1), (left, right));
        assert_eq!(graph.ends(&normalized.graph, 2), (right, antiparticle));
        assert_eq!(
            normalized
                .graph
                .edges()
                .map(|edge| edge.particle)
                .collect::<Vec<_>>(),
            vec![particle(&model, 1); 3],
        );
        assert_eq!(normalized.reversed_edges, BTreeSet::from([EdgeId(1)]));
        assert!(normalized.signs.external_ordering_negative);
        assert!(normalized.signs.include_external_ordering);
        assert!(!normalized.signs.include_antifermion_spin_sum);
        assert!(!normalized.signs.antifermion_spin_sum_negative);

        let symmetrized = ResolvedProcess::new(
            &model,
            &Process::new([1_i64, -1], Vec::<i64>::new()),
            GenerationType::Amplitude,
            &GenerationOptions::default()
                .symmetrize_initial(true)
                .symmetrize_external_fermions(true),
        )
        .unwrap()
        .normalize_fermion_flows(&graph.build(), &model)
        .unwrap();
        assert!(!symmetrized.signs.include_external_ordering);
        assert!(!symmetrized.signs.external_ordering_negative);
    }

    #[test]
    fn ghost_chains_follow_external_anticommuting_order() {
        let source = fermion_model().to_json().unwrap();
        let mut definition = serde_json::from_str::<serde_json::Value>(&source).unwrap();
        for particle in definition["particles"].as_array_mut().unwrap() {
            match particle["name"].as_str() {
                Some("f") => {
                    particle["spin"] = serde_json::json!(-1);
                    particle["ghost_number"] = serde_json::json!(1);
                }
                Some("f~") => {
                    particle["spin"] = serde_json::json!(-1);
                    particle["ghost_number"] = serde_json::json!(-1);
                }
                _ => {}
            }
        }
        definition["lorentz_structures"][0]["spins"] = serde_json::json!([-1, -1, 3]);
        let model = Model::from_json(&definition.to_string()).unwrap();
        let process = ResolvedProcess::new(
            &model,
            &Process::new([1_i64, -1], Vec::<i64>::new()),
            GenerationType::Amplitude,
            &GenerationOptions::default(),
        )
        .unwrap();
        let mut graph = Graph::new();
        let ghost = graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Incoming,
            particle: particle(&model, 1),
        }));
        let left = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let right = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let antighost = graph.add_node(ColoredNode::External(ExternalNode {
            index: 1,
            state: ExternalState::Incoming,
            particle: particle(&model, -1),
        }));
        let flow = particle(&model, 1);
        graph.add_edge(ghost, left, true, flow).unwrap();
        graph.add_edge(right, left, true, flow).unwrap();
        graph.add_edge(right, antighost, true, flow).unwrap();

        let normalized = process
            .normalize_fermion_flows(&graph.build(), &model)
            .unwrap();
        assert_eq!(graph.ends(&normalized.graph, 0), (ghost, left));
        assert_eq!(graph.ends(&normalized.graph, 1), (left, right));
        assert_eq!(graph.ends(&normalized.graph, 2), (right, antighost));
        // The all-incoming order [ghost, antighost] is one transposition
        // from the canonical [antighost, ghost] order, just as for spinors.
        assert!(normalized.signs.external_ordering_negative);
        assert!(!normalized.signs.antifermion_spin_sum_negative);
    }

    /// The `scalars_2p_3p` fixture with the quartic vertex `name` of `particles`.
    fn scalars_with_quartic_vertex(name: &str, particles: [&str; 4]) -> Model {
        let model = Model::from_json(include_str!(
            "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
        ))
        .unwrap();
        let mut definition = serde_json::to_value(&model).unwrap();
        definition["lorentz_structures"]
            .as_array_mut()
            .unwrap()
            .push(serde_json::json!({
                "name":"SCALAR_4_LORENTZ_STRUCTURE", "spins":[1,1,1,1], "structure":"1"
            }));
        definition["vertex_rules"]
            .as_array_mut()
            .unwrap()
            .push(serde_json::json!({
                "name":name, "particles":particles,
                "color_structures":["1"], "lorentz_structures":["SCALAR_4_LORENTZ_STRUCTURE"],
                "couplings":[["SCALAR_COUPLING"]]
            }));
        serde_json::from_value(definition).unwrap()
    }

    #[test]
    fn cp_symmetrization_groups_mirrored_real_scalar_forward_graphs() {
        let model = scalars_with_quartic_vertex("V_4_SCALAR_0000", ["scalar_0"; 4]);
        let mut graphs = Vec::new();
        for mirrored in [false, true] {
            let mut graph = Graph::new();
            // Sew a cubic box to a quartic Born vertex. Its two external
            // neighbors have different valences, so only a global side swap
            // can identify the mirror.
            for node in 0..5 {
                graph.add_node(ColoredNode::Interaction(vertex(
                    &model,
                    if node == 4 {
                        "V_4_SCALAR_0000"
                    } else {
                        "V_3_SCALAR_000"
                    },
                )));
            }
            for index in 0..2 {
                let index = if mirrored { 1 - index } else { index };
                graph.add_node(ColoredNode::External(ExternalNode {
                    index,
                    state: if index == 0 {
                        ExternalState::Incoming
                    } else {
                        ExternalState::Outgoing
                    },
                    particle: particle(&model, 1000),
                }));
            }
            for (left, right) in [
                (5, 0),
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 0),
                (1, 4),
                (2, 4),
                (3, 4),
                (4, 6),
            ] {
                let (left, right) = if mirrored {
                    (right, left)
                } else {
                    (left, right)
                };
                graph
                    .add_edge(left, right, true, edge(&model, 1000))
                    .unwrap();
            }
            graphs.push(graph);
        }
        for enabled in [false, true] {
            let process = ResolvedProcess::new(
                &model,
                &Process::new([1000_i64], [1000_i64, 1000, 1000]),
                GenerationType::CrossSection,
                &GenerationOptions::default()
                    .with_loop_count(3, 3)
                    .unwrap()
                    .symmetrize_left_right(enabled),
            )
            .unwrap();
            let canonical = graphs
                .iter()
                .map(|graph| {
                    process
                        .canonicalize_numerator_graph(
                            &graph.build(),
                            &model,
                            &NumeratorGrouping::None,
                        )
                        .unwrap()
                })
                .collect::<Vec<_>>();
            assert_eq!(canonical[0] == canonical[1], enabled);
        }
    }

    #[test]
    fn global_left_right_canonicalization_preserves_charged_slash_momenta() {
        let model = fermion_model();
        let process = ResolvedProcess::new(
            &model,
            &Process::new([1_i64, -1], [22_i64]),
            GenerationType::CrossSection,
            &GenerationOptions::default().symmetrize_left_right(true),
        )
        .unwrap();
        let mut graph = Graph::new();
        let external_particles = [1_i64, -1, 1, -1];
        let external_nodes = external_particles
            .into_iter()
            .enumerate()
            .map(|(index, pdg)| {
                graph.add_node(ColoredNode::External(ExternalNode {
                    index,
                    state: if index < 2 {
                        ExternalState::Incoming
                    } else {
                        ExternalState::Outgoing
                    },
                    particle: particle(&model, pdg),
                }))
            })
            .collect::<Vec<_>>();
        let left = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let right = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let charged = particle(&model, 1);
        graph
            .add_edge(external_nodes[0], left, true, charged)
            .unwrap();
        graph
            .add_edge(left, external_nodes[1], true, charged)
            .unwrap();
        graph
            .add_edge(right, external_nodes[2], true, charged)
            .unwrap();
        graph
            .add_edge(external_nodes[3], right, true, charged)
            .unwrap();
        graph.add_edge(left, right, true, charged).unwrap();

        let external_swap = BTreeMap::from([
            (external_nodes[0], external_nodes[2]),
            (external_nodes[1], external_nodes[3]),
            (external_nodes[2], external_nodes[0]),
            (external_nodes[3], external_nodes[1]),
        ]);
        let mut mirror = Graph::new();
        for node in &graph.nodes {
            mirror.add_node(*node);
        }
        for &(source, sink, directed, particle) in &graph.edges {
            mirror
                .add_edge(
                    external_swap.get(&source).copied().unwrap_or(source),
                    external_swap.get(&sink).copied().unwrap_or(sink),
                    directed,
                    model.particle_by_id(particle).unwrap().antiparticle,
                )
                .unwrap();
        }

        let grouping = NumeratorGrouping::UpToScalar(GraphGroupingOptions::default());
        let canonical = process
            .canonicalize_numerator_graph(&graph.build(), &model, &grouping)
            .unwrap();
        let canonical_mirror = process
            .canonicalize_numerator_graph(&mirror.build(), &model, &grouping)
            .unwrap();
        assert_eq!(canonical, canonical_mirror);

        let internal_edge = canonical
            .edges()
            .position(|edge| !edge.is_external())
            .unwrap();
        let generator = Generator::new(model.clone());
        let (particle_numerator, particle_fallback) = generator
            .propagator_numerator(model.particle_by_id(particle(&model, 1)).unwrap())
            .unwrap();
        let (antiparticle_numerator, antiparticle_fallback) = generator
            .propagator_numerator(model.particle_by_id(particle(&model, -1)).unwrap())
            .unwrap();
        assert_eq!(particle_numerator, antiparticle_numerator);
        assert!(!particle_fallback);
        assert!(antiparticle_fallback);
        let (propagator, used_antiparticle_fallback) = generator
            .propagator_numerator(
                model
                    .particle_by_id(canonical.edge(internal_edge).particle)
                    .unwrap(),
            )
            .unwrap();
        let (first_flow, second_flow) = if used_antiparticle_fallback {
            (Flow::Sink, Flow::Source)
        } else {
            (Flow::Source, Flow::Sink)
        };
        let legs = [
            NumeratorHalfEdge {
                edge: internal_edge,
                flow: first_flow,
                hedge: 2 * internal_edge + usize::from(first_flow == Flow::Sink),
                spin: 2,
                color: 1,
            },
            NumeratorHalfEdge {
                edge: internal_edge,
                flow: second_flow,
                hedge: 2 * internal_edge + usize::from(second_flow == Flow::Sink),
                spin: 2,
                color: 1,
            },
        ];
        let localized = NumeratorInstantiation::new(NumeratorOwner::Edge(internal_edge), &legs)
            .instantiate(&propagator, NumeratorSector::Spin)
            .unwrap();
        assert_eq!(function_count(&localized, symbol!("gammalooprs::Q")), 1);
    }

    #[test]
    fn canonical_parallel_edges_are_ordered_by_signed_pdg() {
        let model = standard_model();
        let process = ResolvedProcess::new(
            &model,
            &photon_process(),
            GenerationType::Amplitude,
            &GenerationOptions::default(),
        )
        .unwrap();
        let grouping = NumeratorGrouping::UpToSign(GraphGroupingOptions::default());

        for particles in [[-4_i64, 251_i64], [4_i64, -251_i64]] {
            let mut graph = Graph::new();
            let source = graph.add_node(ColoredNode::Interaction(vertex(&model, "V_71")));
            let sink = graph.add_node(ColoredNode::Interaction(vertex(&model, "V_71")));
            for pdg in particles {
                graph
                    .add_edge(source, sink, true, edge(&model, pdg))
                    .unwrap();
            }

            let canonical = process
                .canonicalize_numerator_graph(&graph.build(), &model, &grouping)
                .unwrap();
            let pdgs = canonical
                .edges()
                .map(|edge| model.particle_by_id(edge.particle).unwrap().pdg_code)
                .collect::<Vec<_>>();
            assert!(pdgs.windows(2).all(|pair| pair[0] <= pair[1]), "{pdgs:?}");
        }
    }

    #[test]
    fn cp_transformed_edges_select_the_unique_conjugate_interaction_rule() {
        let model = standard_model();
        let mut graph = Graph::new();
        let interaction = graph.add_node(ColoredNode::Interaction(vertex(&model, "V_82")));
        for (index, pdg) in [-4_i64, 1, 251].into_iter().enumerate() {
            let external = graph.add_node(ColoredNode::External(ExternalNode {
                index,
                state: ExternalState::Incoming,
                particle: particle(&model, pdg),
            }));
            graph
                .add_edge(interaction, external, true, edge(&model, pdg))
                .unwrap();
        }

        let repaired = ResolvedProcess::repair_cp_interaction_rules(graph.build(), &model).unwrap();
        assert_eq!(repaired.vertices().next().unwrap(), vertex(&model, "V_117"));
        let rule = model.vertex_rule("V_117").unwrap();
        assert_eq!(
            Generator::new(model.clone())
                .interaction_half_edges(&repaired, 0, rule)
                .unwrap()
                .len(),
            rule.particles.len()
        );
    }

    #[test]
    fn cp_repaired_interaction_rules_obey_final_vertex_filters() {
        let model = standard_model();
        let mut graph = Graph::new();
        let interaction = graph.add_node(ColoredNode::Interaction(vertex(&model, "V_82")));
        for (index, pdg) in [-4_i64, 1, 251].into_iter().enumerate() {
            let external = graph.add_node(ColoredNode::External(ExternalNode {
                index,
                state: ExternalState::Incoming,
                particle: particle(&model, pdg),
            }));
            graph
                .add_edge(interaction, external, true, edge(&model, pdg))
                .unwrap();
        }
        let repaired = ResolvedProcess::repair_cp_interaction_rules(graph.build(), &model).unwrap();
        assert_eq!(repaired.vertices().next().unwrap(), vertex(&model, "V_117"));

        let resolved_options = |filter| {
            GenerationOptions::default()
                .with_graph_filter(filter)
                .resolve_selectors(&model)
                .unwrap()
        };
        let allow_original =
            resolved_options(GenerationFilter::VertexAllow(vec![VertexSelector::from(
                "V_82",
            )]));
        let veto_repaired =
            resolved_options(GenerationFilter::VertexVeto(vec![VertexSelector::from(
                "V_117",
            )]));
        let allow_repaired =
            resolved_options(GenerationFilter::VertexAllow(vec![VertexSelector::from(
                "V_117",
            )]));
        let generator = Generator::new(model);

        assert!(
            !generator.passes_vertex_selection(&repaired, &allow_original),
            "a CP repair must not bypass an allow-list"
        );
        assert!(
            !generator.passes_vertex_selection(&repaired, &veto_repaired),
            "a CP repair must not bypass a veto-list"
        );
        assert!(
            generator.passes_vertex_selection(&repaired, &allow_repaired),
            "the finalized CP rule should pass when explicitly allowed"
        );
    }

    #[test]
    fn antiparticle_edges_use_exact_propagator_ufo_slots() {
        let source = fermion_model().to_json().unwrap();
        let mut definition = serde_json::from_str::<serde_json::Value>(&source).unwrap();
        let propagators = definition["propagators"].as_array_mut().unwrap();
        propagators
            .iter_mut()
            .find(|propagator| propagator["name"] == "f_prop")
            .unwrap()["numerator"] = serde_json::json!("PSlash(2,1)");
        propagators.push(serde_json::json!({
            "name": "f~_prop",
            "particle": "f~",
            "numerator": "-PSlash(1,2)",
            "denominator": "P^2"
        }));
        let model = Model::from_json(&serde_json::to_string(&definition).unwrap()).unwrap();
        let generator = Generator::new(model.clone());
        let mut graph = Graph::new();
        let interaction = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let external = graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Incoming,
            particle: particle(&model, 22),
        }));
        graph
            .add_edge(interaction, interaction, true, particle(&model, -1))
            .unwrap();
        graph
            .add_edge(external, interaction, false, edge(&model, 22))
            .unwrap();

        let diagram = generator
            .to_diagram(
                Some("antiparticle-propagator".to_owned()),
                GenerationType::Amplitude,
                ColoredTopology {
                    graph: graph.build(),
                    symmetry: 1,
                    cut_partitions: Vec::new(),
                    topology_threshold_candidates: Vec::new(),
                },
                FermionSigns {
                    internal_loop_negative: true,
                    ..FermionSigns::default()
                },
                1,
                &GenerationOptions::default(),
            )
            .unwrap();
        let (actual_edge, _, actual_data) = diagram
            .edges()
            .find(|(_, _, edge)| edge.external.is_none())
            .unwrap();
        let actual = &actual_data.numerator;
        let (antiparticle_propagator, antiparticle_fallback) = generator
            .propagator_numerator(model.particle_by_id(particle(&model, -1)).unwrap())
            .unwrap();
        let (particle_propagator, particle_fallback) = generator
            .propagator_numerator(model.particle_by_id(particle(&model, 1)).unwrap())
            .unwrap();
        assert!(!antiparticle_fallback);
        assert!(!particle_fallback);
        assert_ne!(particle_propagator, antiparticle_propagator);
        let instantiate = |propagator: &Atom, flows: [Flow; 2]| {
            let legs = flows.map(|flow| NumeratorHalfEdge {
                edge: 0,
                flow,
                hedge: usize::from(flow == Flow::Sink),
                spin: 2,
                color: 1,
            });
            NumeratorInstantiation::new(NumeratorOwner::Edge(0), &legs)
                .instantiate(propagator, NumeratorSector::Spin)
                .unwrap()
        };
        let native = |expression: Atom| {
            let source = diagram
                .half_edge_id(DiagramHalfEdge {
                    edge: actual_edge,
                    endpoint: DiagramEndpoint::Source,
                })
                .unwrap();
            let sink = diagram
                .half_edge_id(DiagramHalfEdge {
                    edge: actual_edge,
                    endpoint: DiagramEndpoint::Target,
                })
                .unwrap();
            expression.replace_map_bottom_up(|term, _, output| {
                let AtomView::Fun(function) = term else {
                    return;
                };
                let head = function.get_symbol();
                let arguments = function.iter().collect::<Vec<_>>();
                let owner = usize::try_from(arguments[0]).ok();
                let index = if head == feynkit_graph::symbols::hedge_index() {
                    match owner {
                        Some(0) => Some(source.0),
                        Some(1) => Some(sink.0),
                        _ => None,
                    }
                } else if (head == feynkit_graph::symbols::edge_index()
                    || head == momentum_symbol())
                    && owner == Some(0)
                {
                    Some(actual_edge.0)
                } else {
                    None
                };
                if let Some(index) = index {
                    **output = FunctionBuilder::new(head)
                        .add_arg(index)
                        .add_args(arguments[1..].iter().copied())
                        .finish();
                }
            })
        };
        let exact = native(instantiate(
            &antiparticle_propagator,
            [Flow::Source, Flow::Sink],
        ));
        let blanket_reversal = native(instantiate(
            &antiparticle_propagator,
            [Flow::Sink, Flow::Source],
        ));
        let particle_fallback = native(instantiate(
            &particle_propagator,
            [Flow::Source, Flow::Sink],
        ));

        assert_eq!(actual, &exact);
        assert_ne!(actual, &blanket_reversal);
        assert_ne!(actual, &particle_fallback);
    }

    #[test]
    fn antiparticle_propagator_falls_back_to_particle_definition() {
        let source = fermion_model().to_json().unwrap();
        let mut definition = serde_json::from_str::<serde_json::Value>(&source).unwrap();
        definition["propagators"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|propagator| propagator["name"] == "f_prop")
            .unwrap()["numerator"] = serde_json::json!("PSlash(2,1)");
        let model = Model::from_json(&serde_json::to_string(&definition).unwrap()).unwrap();
        let generator = Generator::new(model.clone());

        let (antiparticle_numerator, antiparticle_fallback) = generator
            .propagator_numerator(model.particle_by_id(particle(&model, -1)).unwrap())
            .unwrap();
        let (particle_numerator, particle_fallback) = generator
            .propagator_numerator(model.particle_by_id(particle(&model, 1)).unwrap())
            .unwrap();
        assert_eq!(antiparticle_numerator, particle_numerator);
        assert!(antiparticle_fallback);
        assert!(!particle_fallback);
    }

    #[test]
    fn vertex_orders_take_each_order_maximum() {
        let model = fermion_model();
        let rule_orders = model.vertex_rule("V").unwrap().coupling_orders(&model);
        assert_eq!(
            rule_orders,
            BTreeMap::from([("QCD".to_owned(), 2), ("QED".to_owned(), 3)])
        );
    }

    fn concrete_interaction_numerator(
        model: &Model,
        rule: VertexRuleId,
        frame: &NumeratorInstantiation<'_>,
    ) -> Result<Atom, GenerationError> {
        let rule = model.vertex_rule_by_id(rule)?;
        let mut terms = Vec::new();
        for (color_index, color) in rule.color_structures.iter().enumerate() {
            for (lorentz_index, lorentz) in rule.lorentz_structures.iter().enumerate() {
                if let Some(coupling) = rule
                    .couplings
                    .get(color_index)
                    .and_then(|row| row.get(lorentz_index))
                    .and_then(Option::as_ref)
                {
                    let coupling = model.coupling_by_id(*coupling)?;
                    let coupling = Atom::var(symbol!(&format!("UFO::{}", coupling.name)));
                    let lorentz = &model.lorentz_structure_by_id(*lorentz)?.structure;
                    let color = frame.instantiate(color, NumeratorSector::Color)?;
                    let lorentz = frame.instantiate(lorentz, NumeratorSector::Spin)?;
                    terms.push(coupling * color * lorentz);
                }
            }
        }
        Ok(if terms.is_empty() {
            Atom::one()
        } else {
            Atom::add_many(terms)
        })
    }

    fn compile_test_factor(expression: &Atom, sector: NumeratorSector) -> InteractionFactor {
        let legs = tensor_test_legs()
            .into_iter()
            .enumerate()
            .map(|(slot, leg)| NumeratorHalfEdge {
                edge: slot,
                hedge: slot,
                ..leg
            })
            .collect::<Vec<_>>();
        InteractionFactor::compile(
            expression,
            sector,
            &NumeratorInstantiation {
                owner: NumeratorOwner::Vertex(0),
                legs: &legs,
                frame: NumeratorFrame::Interaction,
            },
        )
    }

    #[test]
    fn compiled_interaction_slots_preserve_flow_and_self_loop_ports() {
        let mut legs = tensor_test_legs();
        // Sparse native hedge IDs are independent of edge IDs. Both ports of
        // this loop have one momentum, but remain different tensor indices.
        legs[0].edge = 512;
        legs[1].edge = 512;
        legs[0].hedge = 257;
        legs[1].hedge = 65_536;
        legs[2].edge = 255;
        legs[2].hedge = 256;
        let spin = test_atom("(P(3,1)-P(3,2))*Gamma(3,2,1)+P(3)+P()");
        let color = test_atom("T(3,2,1)*f(3,-1,-2)*f(4,-2,-1)");
        for (expression, sector) in [
            (spin, NumeratorSector::Spin),
            (color, NumeratorSector::Color),
        ] {
            let compiled = compile_test_factor(&expression, sector);
            assert!(matches!(compiled, InteractionFactor::Local(_)));
            for source in [true, false] {
                legs[0].flow = if source { Flow::Source } else { Flow::Sink };
                let frame = NumeratorInstantiation::new(NumeratorOwner::Vertex(65_537), &legs);
                assert_eq!(
                    compiled.instantiate(&frame).unwrap(),
                    frame.instantiate(&expression, sector).unwrap()
                );
            }
        }
        let expression = test_atom("P(3,1)+P(3)");
        let compiled = compile_test_factor(&expression, NumeratorSector::Spin);
        legs[0].flow = Flow::Source;
        let frame = NumeratorInstantiation::new(NumeratorOwner::Vertex(99), &legs);
        assert_eq!(
            compiled.instantiate(&frame).unwrap(),
            Atom::Zero,
            "only the explicitly selected momentum receives the incoming-flow sign"
        );
        let gamma = compile_test_factor(&test_atom("Gamma(3,2,1)"), NumeratorSector::Spin)
            .instantiate(&frame)
            .unwrap()
            .to_plain_string();
        assert!(gamma.contains("hedge(257,1)"));
        assert!(gamma.contains("hedge(65536,1)"));
    }

    #[test]
    fn interaction_templates_defer_concrete_walks_and_literal_owners() {
        let mut legs = tensor_test_legs();
        legs[0].edge = 257;
        legs[1].edge = 512;
        legs[0].hedge = 65_536;
        legs[1].hedge = 256;
        for (expression, sector) in [
            (
                "PSlash(1,2,1)+PSlash(1,2,2)+Metric(dummy(4),dummy(4))",
                NumeratorSector::Spin,
            ),
            (
                "Identity(1,-1)*Identity(-1,-2)*T(3,2,-2)",
                NumeratorSector::Color,
            ),
            (
                "Gamma(3,2,1)*gammalooprs::Q(257,spenso::mink(4,gammalooprs::vertex(7,4)))",
                NumeratorSector::Spin,
            ),
            (
                "gammalooprs::hedge(256,1)*gammalooprs::edge(512,1)",
                NumeratorSector::Spin,
            ),
        ] {
            let expression = test_atom(expression);
            let compiled = compile_test_factor(&expression, sector);
            assert!(matches!(compiled, InteractionFactor::Deferred(_, _)));
            for owner in [7, 256, 65_536] {
                let frame = NumeratorInstantiation::new(NumeratorOwner::Vertex(owner), &legs);
                assert_eq!(
                    compiled.instantiate(&frame).unwrap(),
                    frame.instantiate(&expression, sector).unwrap()
                );
            }
        }
        // A momentum used as a color index can acquire a representation
        // conflict only after the two self-loop carriers become equal.
        let aliased_color = test_atom("T(-P(1,1),P(1,2),-1)");
        let compiled = compile_test_factor(&aliased_color, NumeratorSector::Color);
        assert!(matches!(compiled, InteractionFactor::Deferred(_, _)));
        legs[1].edge = legs[0].edge;
        let frame = NumeratorInstantiation::new(NumeratorOwner::Vertex(257), &legs);
        assert_eq!(
            compiled.instantiate(&frame).unwrap_err().to_string(),
            frame
                .instantiate(&aliased_color, NumeratorSector::Color)
                .unwrap_err()
                .to_string()
        );
        let cancelled_error = test_atom("C(1)*(P(3,1)+P(3))");
        let compiled = compile_test_factor(&cancelled_error, NumeratorSector::Spin);
        assert!(matches!(compiled, InteractionFactor::Deferred(_, _)));
        legs[0].flow = Flow::Source;
        let frame = NumeratorInstantiation::new(NumeratorOwner::Vertex(257), &legs);
        assert_eq!(compiled.instantiate(&frame).unwrap(), Atom::Zero);
        legs[0].flow = Flow::Sink;
        let frame = NumeratorInstantiation::new(NumeratorOwner::Vertex(257), &legs);
        assert_eq!(
            compiled.instantiate(&frame).unwrap_err().to_string(),
            frame
                .instantiate(&cancelled_error, NumeratorSector::Spin)
                .unwrap_err()
                .to_string()
        );
        let malformed = test_atom("C(1)");
        let compiled = compile_test_factor(&malformed, NumeratorSector::Spin);
        let frame = NumeratorInstantiation::new(NumeratorOwner::Vertex(257), &legs);
        assert_eq!(
            compiled.instantiate(&frame).unwrap_err().to_string(),
            frame
                .instantiate(&malformed, NumeratorSector::Spin)
                .unwrap_err()
                .to_string()
        );
    }

    #[test]
    fn compiled_interactions_mix_required_operations_and_ignore_unused_factors() {
        let mut definition: serde_json::Value =
            serde_json::from_str(&fermion_model().to_json().unwrap()).unwrap();
        definition["lorentz_structures"].as_array_mut().unwrap().extend([
            serde_json::json!({"name":"FFSlash","spins":[2,2,3],"structure":"PSlash(2,1,1)+PSlash(2,1,2)"}),
            serde_json::json!({"name":"Unused","spins":[2,2,3],"structure":"C(1)"}),
        ]);
        definition["vertex_rules"][0]["lorentz_structures"] =
            serde_json::json!(["FFV", "FFSlash", "Unused"]);
        definition["vertex_rules"][0]["couplings"] =
            serde_json::json!([["GC_LOW", null, null], [null, "GC_HIGH", null]]);
        let model = Model::from_json(&definition.to_string()).unwrap();
        let generator = Generator::new(model.clone());
        let rule = model.vertex_rule_id("V").unwrap();
        let mut legs = tensor_test_legs();
        let first = NumeratorInstantiation::new(NumeratorOwner::Vertex(0), &legs[..3]);
        let compiled = CompiledInteraction::new(&model, rule, &first).unwrap();
        assert!(matches!(compiled, CompiledInteraction::Mixed(_)));
        for index in [1, 257, 65_536] {
            for (slot, leg) in legs.iter_mut().enumerate() {
                leg.edge = index + slot;
                leg.hedge = index + 11 - slot;
            }
            let frame = NumeratorInstantiation::new(NumeratorOwner::Vertex(index), &legs[..3]);
            assert_eq!(
                generator.interaction_numerator(rule, &frame).unwrap(),
                concrete_interaction_numerator(&model, rule, &frame).unwrap()
            );
        }
        assert_eq!(
            generator.numerators.interaction_rules.read().unwrap().len(),
            1,
            "graph positions and momentum flow do not create another compiled rule"
        );
    }

    #[test]
    fn compiled_standard_model_rules_match_concrete_numerators() {
        let model = standard_model();
        let generator = Generator::new(model.clone());
        let mut compiled_count = 0;
        for (rule_index, rule) in model.vertex_rules().iter().enumerate() {
            let rule_id = model.vertex_rule_id_at(rule_index).unwrap();
            let mut legs = rule
                .particles
                .iter()
                .enumerate()
                .map(|(slot, particle)| {
                    let particle = model.particle_by_id(*particle).unwrap();
                    NumeratorHalfEdge {
                        edge: slot,
                        hedge: slot,
                        flow: Flow::Sink,
                        spin: particle.spin,
                        color: particle.color,
                    }
                })
                .collect::<Vec<_>>();
            let local = NumeratorInstantiation::new(NumeratorOwner::Vertex(0), &legs);
            compiled_count += usize::from(matches!(
                CompiledInteraction::new(&model, rule_id, &local).unwrap(),
                CompiledInteraction::Local(_)
            ));
            for (owner, base, loops) in [(0, 0, false), (257, 255, false), (65_537, 65_536, true)] {
                let count = legs.len();
                for (slot, leg) in legs.iter_mut().enumerate() {
                    leg.edge = base + if loops { slot / 2 } else { slot };
                    leg.hedge = base + count - slot;
                    leg.flow = if slot % 2 == 0 {
                        Flow::Source
                    } else {
                        Flow::Sink
                    };
                }
                let frame = NumeratorInstantiation::new(NumeratorOwner::Vertex(owner), &legs);
                assert_eq!(
                    generator.interaction_numerator(rule_id, &frame).unwrap(),
                    concrete_interaction_numerator(&model, rule_id, &frame).unwrap(),
                    "rule {} in owner {owner}, base {base}, loops {loops}",
                    rule.name
                );
            }
        }
        assert_eq!(
            compiled_count,
            model.vertex_rules().len(),
            "the standard-model rules have no graph-dependent lowering operations"
        );
        assert_eq!(
            generator.numerators.interaction_rules.read().unwrap().len(),
            model.vertex_rules().len()
        );
    }

    #[test]
    fn interaction_numerator_keeps_named_couplings_and_lowers_tensors() {
        let model = fermion_model();
        let generator = Generator::new(model.clone());
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(7), &legs[..3]);
        let numerator = generator
            .interaction_numerator(model.vertex_rule_id("V").unwrap(), &instantiation)
            .unwrap();

        let numerator = numerator.to_plain_string();
        assert!(numerator.contains("UFO::GC_LOW"));
        assert!(numerator.contains("UFO::GC_HIGH"));
        assert!(!numerator.contains("g_low"));
        assert!(!numerator.contains("g_high"));
        assert!(numerator.contains("spenso::gamma"));
        assert!(numerator.contains("spenso::t"));
        assert!(!numerator.contains("Gamma("));
        assert!(!numerator.contains("UFO::T("));
    }

    #[test]
    fn lowers_ufo_lorentz_tensors_to_spenso_structures() {
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(7), &legs);
        let gamma = instantiation
            .instantiate(&test_atom("Gamma(3,2,1)"), NumeratorSector::Spin)
            .unwrap();
        assert_eq!(
            gamma,
            test_atom(
                "spenso::gamma(spenso::bis(4,gammalooprs::hedge(41,1)),spenso::bis(4,gammalooprs::hedge(20,1)),spenso::mink(4,gammalooprs::hedge(60,1)))"
            )
        );

        let metric = instantiation
            .instantiate(&test_atom("Metric(3,4)"), NumeratorSector::Spin)
            .unwrap();
        assert_eq!(
            metric,
            test_atom(
                "spenso::g(spenso::mink(4,gammalooprs::hedge(60,1)),spenso::mink(4,gammalooprs::hedge(81,1)))"
            )
        );

        let identity = instantiation
            .instantiate(&test_atom("Identity(1,2)"), NumeratorSector::Spin)
            .unwrap();
        assert_eq!(
            identity,
            test_atom(
                "spenso::g(spenso::bis(4,gammalooprs::hedge(20,1)),spenso::bis(4,gammalooprs::hedge(41,1)))"
            )
        );
    }

    #[test]
    fn lowers_covariant_propagator_with_compact_momentum_slots() {
        use idenso::{CookMode, CookSettings, IndexTooling, tensor::SymbolicTensor};
        use spenso::structure::abstract_index::AbstractIndex;

        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Edge(30), &legs[2..]);
        let numerator = instantiation
            .instantiate(
                &test_atom(
                    "-1i*(Metric(idx(1,1),idx(1,2))-(1-xi)*P(idx(1,1))*P(idx(1,2))/spenso::dot(P(spenso::mink(4)),P(spenso::mink(4))))",
                ),
                NumeratorSector::Spin,
            )
            .unwrap();
        assert_eq!(
            numerator,
            test_atom(
                "-1i*(spenso::g(spenso::mink(4,gammalooprs::hedge(60,1)),spenso::mink(4,gammalooprs::hedge(81,1)))-(1-xi)*gammalooprs::Q(30,spenso::mink(4,gammalooprs::hedge(60,1)))*gammalooprs::Q(30,spenso::mink(4,gammalooprs::hedge(81,1)))/spenso::dot(gammalooprs::Q(30,spenso::mink(4)),gammalooprs::Q(30,spenso::mink(4))))"
            )
        );
        assert_eq!(
            instantiation.lower_spin(numerator.clone()).unwrap(),
            numerator
        );
        assert_eq!(numerator.list_dangling::<AbstractIndex>().unwrap().len(), 2);
        let metric = instantiation
            .instantiate(&test_atom("Metric(1,2)"), NumeratorSector::Spin)
            .unwrap();
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let source = SymbolicTensor::infer(cooking.cook((numerator * metric).as_view())).unwrap();
        let result = source.contract(Default::default()).unwrap();
        let result = result.to_dots().unwrap();
        let scalar = cooking.uncook(result.expression().as_view());
        assert_eq!(
            (scalar - test_atom("-1i*(3+xi)")).together().cancel(),
            Atom::Zero
        );
    }

    #[test]
    fn lowers_ufo_charge_conjugation_and_rejects_invalid_arity() {
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(7), &legs);
        let charge_conjugation = instantiation
            .instantiate(&test_atom("C(1,2)"), NumeratorSector::Spin)
            .unwrap();
        assert_eq!(
            charge_conjugation,
            test_atom(
                "spenso::charge_conjugation(spenso::bis(4,gammalooprs::hedge(20,1)),spenso::bis(4,gammalooprs::hedge(41,1)))"
            )
        );
        assert_eq!(
            instantiation
                .instantiate(&test_atom("C(2,1)"), NumeratorSector::Spin)
                .unwrap(),
            -charge_conjugation,
        );
        for malformed in ["C()", "C(1)", "C(1,2,3)"] {
            assert!(matches!(
                instantiation.instantiate(&test_atom(malformed), NumeratorSector::Spin),
                Err(GenerationError::InvalidNumeratorTensor { .. })
            ));
        }
        assert!(
            instantiation
                .instantiate(&test_atom("C(1,2)"), NumeratorSector::Color)
                .is_err()
        );
    }

    #[test]
    fn lowers_ufo_color_tensors_with_dual_representations() {
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(7), &legs);
        let tensors = instantiation
            .instantiate(
                &test_atom("Identity(1,2)*T(3,2,1)*f(3,-1,-2)"),
                NumeratorSector::Color,
            )
            .unwrap();
        let expected = test_atom(
            "spenso::g(spenso::dind(spenso::cof(3,gammalooprs::hedge(20,1))),spenso::cof(3,gammalooprs::hedge(41,1)))*spenso::t(spenso::coad(8,gammalooprs::hedge(60,1)),spenso::cof(3,gammalooprs::hedge(41,1)),spenso::dind(spenso::cof(3,gammalooprs::hedge(20,1))))*spenso::f(spenso::coad(8,gammalooprs::hedge(60,1)),spenso::coad(8,gammalooprs::vertex(7,1)),spenso::coad(8,gammalooprs::vertex(7,2)))",
        );
        assert_eq!(tensors, expected);
    }

    #[test]
    fn expands_pslash_with_one_shared_collision_free_lorentz_index() {
        let legs = [
            NumeratorHalfEdge {
                edge: 5,
                flow: Flow::Source,
                hedge: 10,
                spin: 2,
                color: 3,
            },
            NumeratorHalfEdge {
                edge: 5,
                flow: Flow::Sink,
                hedge: 11,
                spin: 2,
                color: 3,
            },
        ];
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Edge(5), &legs);
        let slash = instantiation
            .instantiate(
                &test_atom("PSlash(2,1)*Metric(dummy(4),dummy(4))"),
                NumeratorSector::Spin,
            )
            .unwrap();
        assert_eq!(
            slash,
            test_atom(
                "spenso::gamma(spenso::bis(4,gammalooprs::hedge(11,1)),spenso::bis(4,gammalooprs::hedge(10,1)),spenso::mink(4,gammalooprs::edge(5,5)))*gammalooprs::Q(5,spenso::mink(4,gammalooprs::edge(5,5)))*spenso::g(spenso::mink(4,gammalooprs::edge(5,4)),spenso::mink(4,gammalooprs::edge(5,4)))"
            )
        );
    }

    #[test]
    fn lowers_ufo_symmetric_color_tensor_with_standard_normalization() {
        use idenso::{
            CookMode, CookSettings, IndexTooling, color::ColorSimplifySettings,
            tensor::SymbolicTensor,
        };
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let simplify = |expression: Atom| {
            let source = SymbolicTensor::infer(cooking.cook(expression.as_view())).unwrap();
            let result = source
                .simplify_algebra(&idenso::tensor::AlgebraSettings {
                    color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                    ..Default::default()
                })
                .unwrap()
                .contract(idenso::tensor::ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })
                .unwrap();
            cooking.uncook(result.expression().as_view())
        };

        let mut legs = [0, 1, 2].map(|edge| NumeratorHalfEdge {
            edge,
            flow: Flow::Source,
            hedge: 2 * edge,
            spin: 1,
            color: 8,
        });
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(0), &legs);
        let tensor = instantiation
            .instantiate(&test_atom("d(1,2,3)"), NumeratorSector::Color)
            .unwrap();
        assert_eq!(tensor, tensor.spenso_conj());
        assert_eq!(
            tensor,
            instantiation
                .instantiate(&test_atom("d(3,1,2)"), NumeratorSector::Color)
                .unwrap()
        );
        assert_eq!(simplify(tensor.clone() * &tensor), test_atom("40/3"));
        let third = instantiation.index(legs[2], 1).unwrap();
        let other = test_atom("other_color_index");
        let right = tensor.replace(third.to_pattern()).with(other.clone());
        assert_eq!(
            simplify(tensor * right),
            test_atom("5/3")
                * ETS.metric(
                    ColorRepresentation::Adjoint
                        .representation()
                        .to_symbolic([third]),
                    ColorRepresentation::Adjoint
                        .representation()
                        .to_symbolic([other]),
                )
        );
        assert_eq!(
            simplify(
                instantiation
                    .instantiate(&test_atom("d(-1,-1,3)"), NumeratorSector::Color)
                    .unwrap()
            ),
            Atom::Zero
        );
        legs[0].color = 3;
        let incompatible = NumeratorInstantiation::new(NumeratorOwner::Vertex(0), &legs);
        assert!(
            incompatible
                .instantiate(&test_atom("d(1,2,3)"), NumeratorSector::Color)
                .is_err()
        );
    }

    #[test]
    fn lowers_ufo_color_epsilon_and_dual_contractions() {
        use idenso::{CookMode, CookSettings, IndexTooling, tensor::SymbolicTensor};
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let simplify = |expression: Atom| {
            let source = SymbolicTensor::infer(cooking.cook(expression.as_view())).unwrap();
            let result = source
                .simplify_algebra(&idenso::tensor::AlgebraSettings {
                    epsilon: true,
                    ..Default::default()
                })
                .unwrap();
            cooking.uncook(result.expression().as_view())
        };

        for (color, name, conjugate_name) in
            [(3, "Epsilon", "EpsilonBar"), (-3, "EpsilonBar", "Epsilon")]
        {
            let legs = [0, 1, 2].map(|edge| NumeratorHalfEdge {
                edge,
                flow: Flow::Source,
                hedge: 2 * edge,
                spin: 1,
                color,
            });
            let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(0), &legs);
            let tensor = instantiation
                .instantiate(
                    &test_atom(&format!("{name}(1,2,3)")),
                    NumeratorSector::Color,
                )
                .unwrap();
            let exchanged = instantiation
                .instantiate(
                    &test_atom(&format!("{name}(2,1,3)")),
                    NumeratorSector::Color,
                )
                .unwrap();
            assert_eq!(tensor, -exchanged);
            assert_eq!(tensor.spenso_conj().spenso_conj(), tensor);
            assert_eq!(
                simplify(tensor.clone() * tensor.spenso_conj()),
                Atom::num(6)
            );
            assert!(
                instantiation
                    .instantiate(
                        &test_atom(&format!("{conjugate_name}(1,2,3)")),
                        NumeratorSector::Color
                    )
                    .is_err()
            );
            assert!(
                instantiation
                    .instantiate(&test_atom(&format!("{name}(1,2)")), NumeratorSector::Color)
                    .is_err()
            );
            // Identity chains must transmit orientation to the external port,
            // without treating the two ends of a summed index as the same slot.
            let relayed = instantiation
                .instantiate(
                    &test_atom(&format!("Identity(1,-1)*Identity(-1,-2)*{name}(-2,2,3)")),
                    NumeratorSector::Color,
                )
                .unwrap();
            assert_eq!(relayed, tensor);
            let sum = instantiation
                .instantiate(
                    &test_atom(&format!("{name}(1,2,3)+{name}(1,2,-1)*Identity(-1,3)")),
                    NumeratorSector::Color,
                )
                .unwrap();
            assert_eq!(sum, Atom::num(2) * tensor);
        }
        let closed = NumeratorInstantiation::new(NumeratorOwner::Vertex(3), &[]);
        assert_eq!(
            simplify(
                closed
                    .instantiate(
                        &test_atom("Epsilon(-1,-2,-3)*EpsilonBar(-1,-2,-3)"),
                        NumeratorSector::Color,
                    )
                    .unwrap()
            ),
            Atom::num(6)
        );
        assert!(
            closed
                .instantiate(
                    &test_atom("Epsilon(-1,-2,-3)*f(-1,-2,-3)"),
                    NumeratorSector::Color
                )
                .is_err()
        );
    }

    #[test]
    fn lowers_summed_fundamental_color_indices() {
        use idenso::{
            CookMode, CookSettings, color::ColorSimplifySettings, tensor::SymbolicTensor,
        };

        let legs = [0, 1].map(|edge| NumeratorHalfEdge {
            edge,
            flow: Flow::Source,
            hedge: 2 * edge,
            spin: 1,
            color: 8,
        });
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(2), &legs);
        let trace = instantiation
            .instantiate(&test_atom("T(1,-1,-2)*T(2,-2,-1)"), NumeratorSector::Color)
            .unwrap();
        let expected = instantiation
            .instantiate(&test_atom("Identity(1,2)/2"), NumeratorSector::Color)
            .unwrap();
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let source = SymbolicTensor::infer(cooking.cook(trace.as_view())).unwrap();
        let result = source
            .simplify_algebra(&idenso::tensor::AlgebraSettings {
                color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                ..Default::default()
            })
            .unwrap()
            .contract(idenso::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(cooking.uncook(result.expression().as_view()), expected);
    }

    #[test]
    fn keeps_color_and_spin_identity_sectors_distinct() {
        let model = standard_model();
        let generator = Generator::new(model.clone());
        let legs = [
            NumeratorHalfEdge {
                edge: 1,
                flow: Flow::Source,
                hedge: 2,
                spin: 2,
                color: -3,
            },
            NumeratorHalfEdge {
                edge: 2,
                flow: Flow::Sink,
                hedge: 5,
                spin: 2,
                color: 3,
            },
            NumeratorHalfEdge {
                edge: 3,
                flow: Flow::Source,
                hedge: 6,
                spin: 1,
                color: 1,
            },
        ];
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(9), &legs);
        let numerator = generator
            .interaction_numerator(model.vertex_rule_id("V_78").unwrap(), &instantiation)
            .unwrap();

        let numerator = numerator.to_plain_string();
        assert!(numerator.contains("spenso::bis(4"));
        assert!(numerator.contains("spenso::cof(3"));
        assert!(numerator.contains("spenso::dind"));
        assert_eq!(numerator.matches("spenso::g(").count(), 2);
        assert!(!numerator.contains("UFO::Identity("));
    }

    #[test]
    fn coupling_upper_bounds_prune_vertex_signatures_without_applying_lower_bounds() {
        let model = standard_model();
        let options = GenerationOptions::default().with_graph_filter(
            GenerationFilter::CouplingOrders(BTreeMap::from([
                ("QCD".to_owned(), (6, Some(6))),
                ("QED".to_owned(), (0, Some(0))),
            ])),
        );
        let signatures = VertexSignatures::new(&model, &options).unwrap();
        let rules: BTreeSet<_> = signatures.rules.iter().copied().collect();
        assert!(!rules.is_empty());
        for (index, rule) in model.vertex_rules().iter().enumerate() {
            let orders = rule.coupling_orders(&model);
            assert_eq!(
                rules.contains(&model.vertex_rule_id_at(index).unwrap()),
                orders.get("QED").copied().unwrap_or(0) == 0
                    && orders.get("QCD").copied().unwrap_or(0) <= 6
            );
        }
    }

    #[test]
    fn generated_qcd_loop_numerators_are_spenso_ready() {
        let model = standard_model();
        let options = GenerationOptions::default()
            .with_loop_count(1, 1)
            .unwrap()
            .threads(1)
            .max_vertices(2)
            .with_graph_filter(GenerationFilter::CouplingOrders(BTreeMap::from([
                ("QCD".to_owned(), (2, Some(2))),
                ("QED".to_owned(), (0, Some(0))),
            ])))
            .with_graph_filter(GenerationFilter::ParticleVeto(
                vec![-6, -4, -3, -2, -1, 1, 2, 3, 4, 6]
                    .into_iter()
                    .map(ParticleSelector::Pdg)
                    .collect(),
            ));
        let generated = Generator::new(model)
            .generate(
                &Process::new(["g"], ["g"]),
                &options,
                GenerationType::Amplitude,
            )
            .unwrap();
        assert_eq!(generated.diagrams.len(), 3);
        assert_eq!(generated.report.retained_count, 3);
        assert_eq!(
            generated
                .diagrams
                .iter()
                .map(FeynmanDiagram::canonical_key)
                .collect::<Result<BTreeSet<_>, _>>()
                .unwrap()
                .len(),
            3,
            "normalized runtime representatives must remain unique after canonical comparison grouping",
        );

        let residual_heads = [
            "UFO::Gamma(",
            "UFO::Metric(",
            "UFO::PSlash(",
            "UFO::Identity(",
            "UFO::T(",
            "UFO::f(",
        ];
        let mut combined = String::new();
        for diagram in &generated.diagrams {
            combined.push_str(&diagram.numerator().to_plain_string());
            for (_, vertex) in diagram.vertices() {
                combined.push_str(&vertex.numerator.to_plain_string());
            }
            for (_, _, edge) in diagram.edges() {
                combined.push_str(&edge.numerator.to_plain_string());
            }
        }
        for head in residual_heads {
            assert!(!combined.contains(head), "residual tensor head {head}");
        }
        assert!(combined.contains("spenso::gamma("));
        assert!(combined.contains("spenso::t("));
        assert!(combined.contains("spenso::f("));
        assert!(combined.contains("spenso::g("));
    }

    #[test]
    fn dot_rules_keep_explicit_fragments_and_fill_the_others() {
        let model = Arc::new(scalar_model());
        let dot = DotDiagram::from_dot(
            model,
            r#"digraph bubble {
                ext [style=invis];
                a [num="7"];
                ext -> a [particle="phi"];
                a -> b [particle="phi", lmb_id=0, num="11"];
                a -> b [particle="phi"];
                b -> ext [particle="phi"];
            }"#,
        )
        .unwrap();
        let diagram = dot.apply_feynman_rules().unwrap();
        assert_eq!(diagram.numerator(), &(Atom::num(77) * test_atom("GC1")));
        assert_eq!(
            diagram.loop_momentum_basis(),
            dot.diagram.loop_momentum_basis()
        );
        assert_eq!(diagram.projector(), &Atom::one());
        diagram.validate().unwrap();
    }

    /// The external states of `diagram`'s projector as `(head, edge, half-edge)`.
    fn projector_states(diagram: &FeynmanDiagram) -> BTreeSet<(String, usize, usize)> {
        let mut states = BTreeSet::new();
        diagram.projector().as_view().visitor(&mut |term| {
            if let AtomView::Fun(function) = term
                && let [edge, index] = function.iter().collect::<Vec<_>>()[..]
                && let Ok(edge) = usize::try_from(edge)
            {
                let mut hedge = None;
                index.visitor(&mut |inner| {
                    if let AtomView::Fun(inner) = inner
                        && inner.get_symbol() == feynkit_graph::symbols::hedge_index()
                    {
                        hedge = usize::try_from(inner.iter().next().unwrap()).ok();
                    }
                    true
                });
                states.insert((
                    function.get_symbol().get_stripped_name().to_owned(),
                    edge,
                    hedge.unwrap(),
                ));
                return false;
            }
            true
        });
        states
    }

    #[test]
    fn dot_rules_build_an_omitted_projector_in_the_finalized_frame() {
        let model = Arc::new(standard_model());
        let attached = |diagram: &FeynmanDiagram, edge: usize, endpoint| {
            diagram
                .half_edge_id(DiagramHalfEdge {
                    edge: EdgeId(edge),
                    endpoint,
                })
                .unwrap()
                .0
        };
        let amplitude = r#"digraph ee_mumu {
            ext [style=invis];
            ext -> a [id=0, particle="e-"];
            ext -> a [id=1, particle="e+"];
            a -> b [id=2, particle="a"];
            b -> ext [id=3, particle="mu-"];
            b -> ext [id=4, particle="mu+"];
        }"#;
        let diagram = DotDiagram::from_dot(model.clone(), amplitude)
            .unwrap()
            .apply_feynman_rules()
            .unwrap();
        use DiagramEndpoint::{Source, Target};
        assert_eq!(
            projector_states(&diagram),
            BTreeSet::from([
                ("u".to_owned(), 0, attached(&diagram, 0, Target)),
                ("vbar".to_owned(), 1, attached(&diagram, 1, Target)),
                ("ubar".to_owned(), 3, attached(&diagram, 3, Source)),
                ("v".to_owned(), 4, attached(&diagram, 4, Source)),
            ])
        );
        let explicit = DotDiagram::from_dot(
            model.clone(),
            &amplitude.replace(
                "ext [style=invis];",
                r#"graph [projector="2"]; ext [style=invis];"#,
            ),
        )
        .unwrap()
        .apply_feynman_rules()
        .unwrap();
        assert_eq!(explicit.projector(), &Atom::num(2));

        // A sewn initial state carries its incoming copy at its target and its
        // outgoing copy at its source.
        let cross_section = DotDiagram::from_dot(
            model,
            r#"digraph ee_mumu_ee {
                ext [style=invis];
                ext -> a [particle="e-", is_cut=0];
                ext -> a [particle="e+", is_cut=1];
                a -> b [particle="a"];
                b -> c [particle="mu-"];
                b -> c [particle="mu+"];
                c -> d [particle="a"];
                d -> ext [particle="e-", is_cut=0];
                d -> ext [particle="e+", is_cut=1];
            }"#,
        )
        .unwrap()
        .apply_feynman_rules()
        .unwrap();
        assert_eq!(
            projector_states(&cross_section),
            BTreeSet::from([
                ("u".to_owned(), 0, attached(&cross_section, 0, Target)),
                ("ubar".to_owned(), 0, attached(&cross_section, 0, Source)),
                ("vbar".to_owned(), 1, attached(&cross_section, 1, Target)),
                ("v".to_owned(), 1, attached(&cross_section, 1, Source)),
            ])
        );
        cross_section.validate().unwrap();
    }

    #[test]
    fn apply_feynman_rules_to_compact_dot() {
        let generator = Generator::new(scalar_model());
        let diagram = FeynmanDiagram::from_dot(
            Arc::clone(&generator.model),
            r#"digraph bubble {
                ext [style=invis];
                ext -> a [particle="phi"];
                a -> b [particle="phi", lmb_id=0];
                a -> b [particle="phi"];
                b -> ext [particle="phi"];
            }"#,
        )
        .unwrap();
        assert_eq!(diagram.numerator(), &Atom::one());
        let routed = diagram.loop_momentum_basis().clone();
        let id = diagram.id().unwrap();
        let diagram = diagram.apply_feynman_rules().unwrap();
        assert_eq!(diagram.numerator(), &test_atom("GC1^2"));
        assert_eq!(diagram.loop_momentum_basis(), &routed);
        assert_eq!(diagram.id().unwrap(), id);
        diagram.validate().unwrap();

        let invalid_slots = diagram
            .map_data(
                |_, vertex| vertex.clone(),
                |_, _, edge| edge.clone().with_slots(VertexSlot(7), VertexSlot(7)),
            )
            .unwrap();
        assert!(matches!(
            invalid_slots.apply_feynman_rules(),
            Err(GenerationError::Diagram(
                DiagramError::DuplicateVertexSlot { .. }
            ))
        ));
    }

    #[test]
    fn numerator_binding_uses_ufo_slots_and_both_self_loop_ends() {
        let generator = Generator::new(standard_model());
        let rule = vertex(&generator.model, "V_98"); // e+, e-, a
        let electron = particle(&generator.model, 11);
        let photon = particle(&generator.model, 22);
        let mut builder = FeynmanDiagramBuilder::new(Arc::clone(&generator.model), None);
        let interaction = builder.add_vertex(DiagramVertex::interaction("v", rule));
        let spectator = builder.add_vertex(DiagramVertex {
            name: "spectator".into(),
            interaction: None,
            numerator: Atom::num(5),
        });

        // Insert slot 2 first: incidence order is not UFO leg order.
        let mut external = DiagramEdge::new(photon, false);
        external.external = Some(ExternalLeg {
            name: "q0".into(),
            index: 0,
            state: ExternalState::Incoming,
            connection: 0,
        });
        builder
            .add_edge_with_slots(None, interaction, external, VertexSlot(0), VertexSlot(2))
            .unwrap();
        builder
            .add_edge_with_slots(
                interaction,
                interaction,
                DiagramEdge::new(electron, true),
                VertexSlot(0),
                VertexSlot(1),
            )
            .unwrap();
        // A dummy occupies a slot beyond the interaction legs and must not bind one.
        let mut dummy = DiagramEdge::new(photon, false);
        dummy.is_dummy = true;
        dummy.numerator = Atom::num(7);
        builder
            .add_edge_with_slots(interaction, spectator, dummy, VertexSlot(3), VertexSlot(0))
            .unwrap();
        let diagram = builder.build().unwrap();
        let (interaction, _) = diagram
            .vertices()
            .find(|(_, data)| data.interaction.is_some())
            .unwrap();
        let (external, _, _) = diagram
            .edges()
            .find(|(_, _, data)| data.external.is_some())
            .unwrap();
        let (internal, _, _) = diagram
            .edges()
            .find(|(_, _, data)| data.external.is_none() && !data.is_dummy)
            .unwrap();
        let legs = [
            (internal, DiagramEndpoint::Source, Flow::Source, 2),
            (internal, DiagramEndpoint::Target, Flow::Sink, 2),
            (external, DiagramEndpoint::Target, Flow::Sink, 3),
        ]
        .map(|(edge, endpoint, flow, spin)| NumeratorHalfEdge {
            edge: edge.0,
            hedge: diagram
                .half_edge_id(DiagramHalfEdge { edge, endpoint })
                .unwrap()
                .0,
            flow,
            spin,
            color: 1,
        });
        assert_ne!(legs[0].hedge, legs[1].hedge);
        let expected = generator
            .interaction_numerator(
                rule,
                &NumeratorInstantiation::new(NumeratorOwner::Vertex(interaction.0), &legs),
            )
            .unwrap();
        let diagram = diagram.apply_feynman_rules().unwrap();
        assert_eq!(diagram.vertex(interaction).unwrap().numerator, expected);
        assert_eq!(
            diagram
                .vertices()
                .find(|(_, data)| data.interaction.is_none())
                .unwrap()
                .1
                .numerator,
            Atom::num(5)
        );
        assert_eq!(
            diagram
                .edges()
                .find(|(_, _, data)| data.is_dummy)
                .unwrap()
                .2
                .numerator,
            Atom::num(7)
        );
    }

    #[test]
    fn apply_feynman_rules_matches_generated_local_factors() {
        let generator = Generator::new(standard_model());
        let generated = generator
            .generate(
                &Process::new(["g"], ["g"]),
                &GenerationOptions::default()
                    .with_loop_count(1, 1)
                    .unwrap()
                    .threads(1)
                    .max_vertices(2)
                    .with_graph_filter(GenerationFilter::CouplingOrders(BTreeMap::from([
                        ("QCD".to_owned(), (2, Some(2))),
                        ("QED".to_owned(), (0, Some(0))),
                    ]))),
                GenerationType::Amplitude,
            )
            .unwrap();
        assert!(!generated.diagrams.is_empty());
        for diagram in generated.diagrams {
            let topology = diagram.clone().with_numerator(Atom::one()).unwrap();
            let restored = topology.apply_feynman_rules().unwrap();
            assert_eq!(restored.to_json().unwrap(), diagram.to_json().unwrap());
            let reapplied = restored.apply_feynman_rules().unwrap();
            assert_eq!(reapplied.to_json().unwrap(), restored.to_json().unwrap());
        }
    }

    #[test]
    fn closed_fermion_loops_always_contribute_their_sign() {
        let generator = Generator::new(fermion_model());
        let mut graph = Graph::new();
        let interaction = graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V")));
        let external = graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Outgoing,
            particle: particle(&generator.model, 22),
        }));
        graph
            .add_edge(
                interaction,
                interaction,
                true,
                particle(&generator.model, 1),
            )
            .unwrap();
        graph
            .add_edge(interaction, external, false, particle(&generator.model, 22))
            .unwrap();

        let loop_count = generator
            .closed_fermion_loop_count(&graph.build(), false)
            .unwrap();
        assert_eq!(loop_count, 1);
        let diagram = generator
            .to_diagram(
                Some("loop".to_owned()),
                GenerationType::Amplitude,
                ColoredTopology {
                    graph: graph.build(),
                    symmetry: 3,
                    cut_partitions: Vec::new(),
                    topology_threshold_candidates: Vec::new(),
                },
                FermionSigns {
                    internal_loop_negative: loop_count % 2 == 1,
                    ..FermionSigns::default()
                },
                1,
                &GenerationOptions::default(),
            )
            .unwrap();
        assert_eq!(
            diagram.overall_factor(),
            &factor_atom("InternalFermionLoopSign(-1)/AutG(3)")
        );
    }

    #[test]
    fn open_ghost_chain_exchange_has_grassmann_sign() {
        let model = standard_model();
        let process = ResolvedProcess::new(
            &model,
            &Process::new([9000005_i64, 9000005], [9000005_i64, 9000005]),
            GenerationType::Amplitude,
            &GenerationOptions::default(),
        )
        .unwrap();
        let mut signs = Vec::new();
        for outgoing in [[2, 3], [3, 2]] {
            let mut graph = Graph::new();
            for index in 0..4 {
                graph.add_node(ColoredNode::External(ExternalNode {
                    index,
                    state: if index < 2 {
                        ExternalState::Incoming
                    } else {
                        ExternalState::Outgoing
                    },
                    particle: particle(&model, 9000005),
                }));
            }
            for _ in 0..2 {
                graph.add_node(ColoredNode::Interaction(vertex(&model, "V_35")));
            }
            for (left, right, pdg) in [
                (0, 4, 9000005),
                (4, outgoing[0], 9000005),
                (1, 5, 9000005),
                (5, outgoing[1], 9000005),
                (4, 5, 21),
            ] {
                graph
                    .add_edge(left, right, true, particle(&model, pdg))
                    .unwrap();
            }
            signs.push(
                process
                    .normalize_fermion_flows(&graph.build(), &model)
                    .unwrap()
                    .signs
                    .external_ordering_negative,
            );
        }
        assert_ne!(
            signs[0], signs[1],
            "exchanging two identical external ghosts must reverse the sign"
        );
    }

    #[test]
    fn closed_ghost_loops_contribute_a_grassmann_minus_sign() {
        let generator = Generator::new(standard_model());
        let generated = generator
            .generate(
                &Process::new(["g"], ["g"]),
                &GenerationOptions::default()
                    .with_loop_count(1, 1)
                    .unwrap()
                    .threads(1)
                    .max_vertices(2)
                    .with_graph_filter(GenerationFilter::VertexAllow(vec![VertexSelector::Name(
                        "V_35".to_owned(),
                    )])),
                GenerationType::Amplitude,
            )
            .unwrap();
        assert_eq!(generated.diagrams.len(), 1);
        let diagram = &generated.diagrams[0];
        assert_eq!(
            diagram
                .edges()
                .filter(|(_, _, edge)| {
                    generator.model.particle_by_id(edge.particle).unwrap().spin == -1
                })
                .count(),
            2,
        );

        // The ghost-gluon vertex and propagator supply no closed-loop sign.
        let factor = diagram.overall_factor();
        let without_loop_sign = factor
            .replace(factor_atom("InternalFermionLoopSign(-1)").to_pattern())
            .with(Atom::one().to_pattern());
        assert_ne!(
            factor, &without_loop_sign,
            "a closed ghost loop must carry its Grassmann minus sign in the overall factor",
        );
    }

    /// Known pre-existing failure: building the half-integer-ghost fixture fails model
    /// validation ("Covariant cut multiplet -24 must contain one massive vector, its
    /// Goldstone, ghost and antighost with the same symbolic pole mass, charge and color"),
    /// before any generation code runs. `failing::` tests stay out of the curated suites.
    mod failing {
        use super::*;

        #[test]
        fn closed_anticommutating_loop_counts() {
            let model = standard_model();
            let mut definition = serde_json::to_value(&model).unwrap();
            for particle in definition["particles"].as_array_mut().unwrap() {
                if particle["ghost_number"].as_i64().unwrap_or(0) != 0 {
                    particle["spin"] = serde_json::json!(2);
                }
            }
            for lorentz in definition["lorentz_structures"].as_array_mut().unwrap() {
                for spin in lorentz["spins"].as_array_mut().unwrap() {
                    if spin.as_i64() == Some(-1) {
                        *spin = serde_json::json!(2);
                    }
                }
            }
            let half_integer_ghost_model: Model = serde_json::from_value(definition).unwrap();
            let generators = [
                Generator::new(model),
                Generator::new(half_integer_ghost_model),
            ];
            for (name, edges, expected, expected_fermions) in [
                ("ghost", vec![(0, 1, 9000005), (1, 0, 9000005)], 1, 0),
                ("antighost", vec![(0, 1, -9000005), (1, 0, -9000005)], 1, 0),
                ("ghost self-loop", vec![(0, 0, 9000005)], 1, 0),
                ("Dirac", vec![(0, 1, 5), (1, 0, 5)], 1, 1),
                ("anti-Dirac", vec![(0, 1, -5), (1, 0, -5)], 1, 1),
                ("Dirac self-loop", vec![(0, 0, 5)], 1, 1),
                ("open Dirac chain", vec![(0, 1, 5), (1, 2, 5)], 0, 0),
                (
                    "open ghost chain",
                    vec![(0, 1, 9000005), (1, 2, 9000005)],
                    0,
                    0,
                ),
                ("gluon cycle", vec![(0, 1, 21), (1, 0, 21)], 0, 0),
                (
                    "ghost and Dirac loops joined by a gluon",
                    vec![
                        (0, 1, 9000005),
                        (1, 0, 9000005),
                        (1, 2, 21),
                        (2, 3, 5),
                        (3, 2, 5),
                    ],
                    2,
                    1,
                ),
                (
                    "two ghost loops joined by a gluon",
                    vec![
                        (0, 1, 9000005),
                        (1, 0, 9000005),
                        (1, 2, 21),
                        (2, 3, 9000005),
                        (3, 2, 9000005),
                    ],
                    2,
                    0,
                ),
                (
                    "two Dirac loops joined by a gluon",
                    vec![(0, 1, 5), (1, 0, 5), (1, 2, 21), (2, 3, 5), (3, 2, 5)],
                    2,
                    2,
                ),
            ] {
                for generator in &generators {
                    let model = &generator.model;
                    let mut graph = Graph::new();
                    let last_node = edges.iter().flat_map(|(a, b, _)| [*a, *b]).max().unwrap();
                    for _ in 0..=last_node {
                        graph.add_node(ColoredNode::Interaction(vertex(model, "V_35")));
                    }
                    for &(a, b, pdg) in &edges {
                        graph.add_edge(a, b, true, edge(model, pdg)).unwrap();
                    }
                    assert_eq!(
                        generator
                            .closed_fermion_loop_count(&graph.build(), true)
                            .unwrap(),
                        expected,
                        "{name}"
                    );
                    assert_eq!(
                        generator
                            .closed_fermion_loop_count(&graph.build(), false)
                            .unwrap(),
                        expected_fermions,
                        "fermion filter: {name}"
                    );
                }
            }
        }
    }

    #[test]
    fn ghost_loop_count_excludes_open_chains_and_counts_each_closed_component() {
        let generator = Generator::new(standard_model());
        let mut graph = Graph::new();
        let incoming = graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Incoming,
            particle: particle(&generator.model, 9000005),
        }));
        let outgoing = graph.add_node(ColoredNode::External(ExternalNode {
            index: 1,
            state: ExternalState::Outgoing,
            particle: particle(&generator.model, 9000005),
        }));
        let gluon = graph.add_node(ColoredNode::External(ExternalNode {
            index: 2,
            state: ExternalState::Outgoing,
            particle: particle(&generator.model, 21),
        }));
        // an open ghost line through one vertex
        let line = graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V_35")));
        let ghost = edge(&generator.model, 9000005);
        graph.add_edge(incoming, line, true, ghost).unwrap();
        graph.add_edge(line, outgoing, true, ghost).unwrap();
        graph
            .add_edge(line, gluon, false, edge(&generator.model, 21))
            .unwrap();
        assert_eq!(
            generator
                .closed_fermion_loop_count(&graph.build(), true)
                .unwrap(),
            0
        );

        for loop_count in 1..=2 {
            let interaction =
                graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V_35")));
            graph
                .add_edge(
                    interaction,
                    interaction,
                    true,
                    edge(&generator.model, 9000005),
                )
                .unwrap();
            assert_eq!(
                generator
                    .closed_fermion_loop_count(&graph.build(), true)
                    .unwrap(),
                loop_count,
            );
            assert_eq!(
                generator
                    .closed_fermion_loop_count(&graph.build(), false)
                    .unwrap(),
                0
            );
        }
    }

    #[test]
    fn fermion_loop_filter_excludes_ghost_loops() {
        let generator = Generator::new(standard_model());
        for (fermion_loops, expected_diagrams) in [(0, 1), (1, 0)] {
            let generated = generator
                .generate(
                    &Process::new(["g"], ["g"]),
                    &GenerationOptions::default()
                        .with_loop_count(1, 1)
                        .unwrap()
                        .threads(1)
                        .max_vertices(2)
                        .with_graph_filter(GenerationFilter::VertexAllow(vec![
                            VertexSelector::Name("V_35".to_owned()),
                        ]))
                        .with_graph_filter(GenerationFilter::FermionLoopCountRange((
                            fermion_loops,
                            fermion_loops,
                        ))),
                    GenerationType::Amplitude,
                )
                .unwrap();
            assert_eq!(generated.diagrams.len(), expected_diagrams);
        }
    }

    #[test]
    fn positive_external_fermion_signs_keep_symbolic_provenance() {
        let generator = Generator::new(fermion_model());
        let diagram = generator
            .to_diagram(
                Some("positive".to_owned()),
                GenerationType::Amplitude,
                ColoredTopology {
                    graph: GeneratedGraph::new([], []),
                    symmetry: 1,
                    cut_partitions: Vec::new(),
                    topology_threshold_candidates: Vec::new(),
                },
                FermionSigns {
                    internal_loop_negative: false,
                    include_external_ordering: true,
                    external_ordering_negative: false,
                    include_antifermion_spin_sum: true,
                    antifermion_spin_sum_negative: false,
                },
                1,
                &GenerationOptions::default(),
            )
            .unwrap();
        assert_eq!(
            diagram.overall_factor(),
            &factor_atom("AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*1/AutG(1)")
        );
    }

    #[test]
    fn branching_fermion_chains_are_rejected() {
        let generator = Generator::new(fermion_model());
        let mut graph = Graph::new();
        let center = graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V")));
        for _ in 0..3 {
            let endpoint = graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V")));
            graph
                .add_edge(center, endpoint, true, particle(&generator.model, 1))
                .unwrap();
        }

        assert!(matches!(
            generator.closed_fermion_loop_count(&graph.build(), false),
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
            generator.generate(&process, &unsupported, GenerationType::Amplitude),
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
            generator.generate(&process, &perturbative, GenerationType::Amplitude),
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
                &GenerationOptions::default().with_graph_filter(sewn_filter),
                GenerationType::Amplitude
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
            generator.generate(&process, &duplicate, GenerationType::Amplitude),
            Err(GenerationError::DuplicateFilter {
                scope: FilterScope::Graph,
                filter: GenerationFilterKind::LoopCountRange
            })
        ));

        let invalid = GenerationOptions::default()
            .with_graph_filter(GenerationFilter::LoopCountRange((2, 1)));
        assert!(matches!(
            generator.generate(&process, &invalid, GenerationType::Amplitude),
            Err(GenerationError::InvalidFilterRange {
                scope: FilterScope::Graph,
                filter: GenerationFilterKind::LoopCountRange,
                minimum: 2,
                maximum: 1
            })
        ));

        let invalid_scope = GenerationOptions::default().with_cut_amplitude_filter(
            GenerationFilter::ParticleVeto(vec![ParticleSelector::Pdg(1)]),
        );
        assert!(matches!(
            generator.generate(&process, &invalid_scope, GenerationType::Amplitude),
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

        let grouping = GenerationOptions::default().filter_zero_color(true);
        generator
            .validate_options(GenerationType::Amplitude, &grouping)
            .unwrap();

        let invalid = GenerationOptions::default().numerator_grouping(
            NumeratorGrouping::Identical(GraphGroupingOptions {
                number_of_numerical_samples: 0,
                ..GraphGroupingOptions::default()
            }),
        );
        assert!(matches!(
            generator.validate_options(GenerationType::Amplitude, &invalid),
            Err(GenerationError::InvalidNumericalSampleCount)
        ));
    }

    #[test]
    fn public_limits_report_arithmetic_overflow() {
        let model = fermion_model();
        let options =
            GenerationOptions::default().with_graph_filter(GenerationFilter::PerturbativeOrders(
                BTreeMap::from([("QED".to_owned(), usize::MAX), ("QCD".to_owned(), 1)]),
            ));
        assert!(matches!(
            unresolved_cut_content(&model, &options),
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
        let unresolved = unresolved_cut_content(&model, &options).unwrap().unwrap();
        assert_eq!(unresolved.maximum_multiplicity, 1);
        assert_eq!(
            unresolved.particles,
            BTreeSet::from([
                particle(&model, -1),
                particle(&model, 1),
                particle(&model, 22),
            ])
        );
        assert!(ResolvedProcess::matches_cut_content(
            &[particle(&model, 1), particle(&model, 22)],
            &[particle(&model, 1)],
            Some(&unresolved)
        ));
        assert!(!ResolvedProcess::matches_cut_content(
            &[
                particle(&model, 1),
                particle(&model, 22),
                particle(&model, 22),
            ],
            &[particle(&model, 1)],
            Some(&unresolved)
        ));
        assert!(!ResolvedProcess::matches_cut_content(
            &[particle(&model, 1), particle(&model, 25)],
            &[particle(&model, 1)],
            Some(&unresolved)
        ));
    }

    #[test]
    fn exact_cut_rejects_species_only_false_positives() {
        let model = exact_cut_model();
        let interaction = |name: &str| ColoredNode::Interaction(vertex(&model, name));
        let external = |index, state, pdg| {
            ColoredNode::External(ExternalNode {
                index,
                state,
                particle: particle(&model, pdg),
            })
        };
        let edge = |pdg| edge(&model, pdg);

        let mut three_loop = Graph::new();
        let h0 = three_loop.add_node(interaction("V_78"));
        let h1 = three_loop.add_node(interaction("V_78"));
        let g0 = three_loop.add_node(interaction("V_76"));
        let p0 = three_loop.add_node(interaction("V_73"));
        let g1 = three_loop.add_node(interaction("V_76"));
        let p1 = three_loop.add_node(interaction("V_73"));
        let q0 = three_loop.add_node(interaction("V_98"));
        let q1 = three_loop.add_node(interaction("V_98"));
        let x0 = three_loop.add_node(external(0, ExternalState::Incoming, -11));
        let x1 = three_loop.add_node(external(1, ExternalState::Incoming, 11));
        let x2 = three_loop.add_node(external(2, ExternalState::Outgoing, -11));
        let x3 = three_loop.add_node(external(3, ExternalState::Outgoing, 11));
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
        let x0 = four_loop.add_node(external(0, ExternalState::Incoming, -11));
        let x1 = four_loop.add_node(external(1, ExternalState::Incoming, 11));
        let x2 = four_loop.add_node(external(2, ExternalState::Outgoing, -11));
        let x3 = four_loop.add_node(external(3, ExternalState::Outgoing, 11));
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
            &Process::new([-11_i64, 11], [5_i64, -5, 25]),
            GenerationType::CrossSection,
            &GenerationOptions::default(),
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
                    .cut_partitions(&graph.build(), &model, &options)
                    .unwrap()
                    .physical
                    .is_empty()
            );
        }
    }

    #[test]
    fn sewn_filter_vetoes_factorization_created_by_pairwise_sewing() {
        let model = fermion_model();
        let process = ResolvedProcess::new(
            &model,
            &Process::new([22_i64, 22], [22_i64]),
            GenerationType::CrossSection,
            &GenerationOptions::default(),
        )
        .unwrap();
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
                    particle: particle(&model, 22),
                }))
            })
            .collect();
        let center = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let photon = edge(&model, 22);
        for external in &externals {
            graph.add_edge(*external, center, false, photon).unwrap();
        }

        let options = GenerationOptions::default()
            .with_graph_filter(GenerationFilter::Sewn(SewnFilterOptions::default()));
        assert!(process.passes_sewn_filter(&graph.build(), &options));

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
                    particle: particle(&model, 22),
                }))
            })
            .collect();
        let left = factorized.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        let right = factorized.add_node(ColoredNode::Interaction(vertex(&model, "V")));
        for external in [externals[0], externals[2]] {
            factorized.add_edge(external, left, false, photon).unwrap();
        }
        for external in [externals[1], externals[3]] {
            factorized.add_edge(external, right, false, photon).unwrap();
        }
        factorized.add_edge(left, right, false, photon).unwrap();
        assert!(!process.passes_sewn_filter(&factorized.build(), &options));
        assert!(process.passes_sewn_filter(
            &factorized.build(),
            &GenerationOptions::default().with_graph_filter(GenerationFilter::Sewn(
                SewnFilterOptions {
                    filter_tadpoles: false,
                }
            ))
        ));
    }

    /// The cross sections of `g g -> u u~` at one loop: their IDs, the sum of their evaluated
    /// overall factors, and the boxes among them as
    /// `(sewn copies on adjacent vertices, evaluated overall factor)`.
    fn gluon_fusion_boxes(symmetrize_initial: bool) -> (Vec<DiagramId>, Atom, Vec<(bool, Atom)>) {
        let process = Process::new(["g", "g"], ["u", "u~"])
            .with_particle_filters(
                Some(["d", "c", "s", "t"].into_iter().map(Into::into).collect()),
                None,
                None,
                vec![],
            )
            .unwrap();
        let options = GenerationOptions::default()
            .with_loop_count(1, 1)
            .unwrap()
            .threads(1)
            .symmetrize_initial(symmetrize_initial);
        let result = Generator::new(Model::qcd())
            .generate(&process, &options, GenerationType::CrossSection)
            .unwrap();
        let mut boxes = vec![];
        for diagram in &result.diagrams {
            let mut internal = BTreeSet::new();
            let mut sewn = vec![];
            for (_, ends, edge) in diagram.edges() {
                let ends = (ends.source.unwrap(), ends.target.unwrap());
                if edge.external.is_some() {
                    sewn.push(ends);
                } else if edge.directed {
                    internal.insert(ends);
                    internal.insert((ends.1, ends.0));
                }
            }
            // the quark loop through all four vertices
            if diagram.vertices().count() == 4 && internal.len() == 8 {
                let adjacent = sewn
                    .iter()
                    .map(|ends| internal.contains(ends))
                    .collect::<Vec<_>>();
                assert_eq!(adjacent[0], adjacent[1]);
                boxes.push((adjacent[0], evaluated_overall_factor(diagram)));
            }
        }
        boxes.sort();
        (
            result
                .diagrams
                .iter()
                .map(|diagram| diagram.id().unwrap())
                .collect(),
            Atom::add_many(result.diagrams.iter().map(evaluated_overall_factor)),
            boxes,
        )
    }

    fn evaluated_overall_factor(diagram: &FeynmanDiagram) -> Atom {
        feynkit_graph::expressions::evaluate_overall_factor(diagram.overall_factor().as_view())
    }

    #[test]
    fn symmetrized_cross_sections_permute_sewn_pairs_together() {
        let (exact_ids, exact_total, exact_boxes) = gluon_fusion_boxes(false);
        assert_eq!(exact_ids.len(), 9);
        assert_eq!(exact_total, Atom::num(-9));
        // the crossed and the planar sewing, each with both orientations of the quark loop
        let box_factor = |adjacent, factor| (adjacent, Atom::num(factor));
        assert_eq!(
            exact_boxes,
            [
                box_factor(false, -1),
                box_factor(false, -1),
                box_factor(true, -1),
                box_factor(true, -1)
            ]
        );
        let (ids, total, boxes) = gluon_fusion_boxes(true);
        // exchanging both sewn pairs reverses the quark loop: one box of each sewing remains,
        // standing for both orientations
        assert_eq!(boxes, [box_factor(false, -2), box_factor(true, -2)]);
        assert_eq!(ids.len(), 5);
        assert_eq!(total, exact_total);
    }

    /// The sorted evaluated overall factors of a one-loop cross section, the tree-level cuts of
    /// a `2 -> 2` process.
    fn born_cross_section_factors(
        model: &Model,
        process: &Process,
        symmetrize_initial: bool,
    ) -> Result<Vec<Atom>, GenerationError> {
        let options = GenerationOptions::default()
            .with_loop_count(1, 1)
            .unwrap()
            .threads(1)
            .symmetrize_initial(symmetrize_initial);
        let mut factors = Generator::new(model.clone())
            .generate(process, &options, GenerationType::CrossSection)?
            .diagrams
            .iter()
            .map(evaluated_overall_factor)
            .collect::<Vec<_>>();
        factors.sort();
        Ok(factors)
    }

    #[test]
    fn beam_exchange_weights_each_kept_graph_by_its_orbit() {
        let model = scalars_with_quartic_vertex(
            "V_4_SCALAR_1122",
            ["scalar_1", "scalar_1", "scalar_2", "scalar_2"],
        );
        let process = |incoming: [&str; 2], outgoing: [&str; 2], vertices: &[&str]| {
            Process::new(incoming, outgoing)
                .with_particle_filters(
                    None,
                    None,
                    Some(vertices.iter().map(|&vertex| vertex.into()).collect()),
                    vec![],
                )
                .unwrap()
        };
        let factors = |values: &[(i64, i64)]| {
            values
                .iter()
                .map(|&value| Atom::num(value))
                .collect::<Vec<_>>()
        };
        for (process, exact, symmetrized) in [
            // The contact graph is its own beam exchange.
            (
                process(
                    ["scalar_1", "scalar_1"],
                    ["scalar_2", "scalar_2"],
                    &["V_4_SCALAR_1122"],
                ),
                factors(&[(1, 2)]),
                factors(&[(1, 2)]),
            ),
            // So are the squared t-channel and the t/u interference graph.
            (
                process(
                    ["scalar_1", "scalar_1"],
                    ["scalar_2", "scalar_2"],
                    &["V_3_SCALAR_112"],
                ),
                factors(&[(1, 1), (1, 1)]),
                factors(&[(1, 1), (1, 1)]),
            ),
            // Distinct final states: the exchange maps |M_t|^2 to |M_u|^2 and M_t M_u* to
            // M_u M_t*, so one graph of each pair remains with weight two.
            (
                process(
                    ["scalar_1", "scalar_1"],
                    ["scalar_0", "scalar_2"],
                    &["V_3_SCALAR_011", "V_3_SCALAR_112"],
                ),
                factors(&[(1, 1); 4]),
                factors(&[(2, 1); 2]),
            ),
            // Distinct beams have no exchange.
            (
                process(
                    ["scalar_1", "scalar_2"],
                    ["scalar_1", "scalar_2"],
                    &["V_4_SCALAR_1122"],
                ),
                factors(&[(1, 1)]),
                factors(&[(1, 1)]),
            ),
        ] {
            assert_eq!(
                born_cross_section_factors(&model, &process, false).unwrap(),
                exact
            );
            assert_eq!(
                born_cross_section_factors(&model, &process, true).unwrap(),
                symmetrized
            );
        }
    }

    #[test]
    fn beam_exchange_keeps_fermion_signs_and_exact_totals() {
        let factors = |values: &[i64]| {
            values
                .iter()
                .map(|&value| Atom::num(value))
                .collect::<Vec<_>>()
        };
        // Møller scattering through photon exchange: the squared and the interference graph,
        // each its own beam exchange, with the relative Fermi sign.
        let moller = Process::new(["e-", "e-"], ["e-", "e-"]);
        for symmetrize_initial in [false, true] {
            assert_eq!(
                born_cross_section_factors(&Model::qed(), &moller, symmetrize_initial).unwrap(),
                factors(&[-1, 1])
            );
        }
        // a a -> t t~: |M_t|^2 and |M_u|^2 form one exchange orbit, M_t M_u* and M_u M_t* the
        // other.
        let photon_fusion = Process::new(["a", "a"], ["t", "t~"])
            .with_particle_filters(None, Some(vec!["a".into(), "t".into()]), None, vec![])
            .unwrap();
        let model = Model::standard_model();
        assert_eq!(
            born_cross_section_factors(&model, &photon_fusion, false).unwrap(),
            factors(&[-1, -1, -1, -1])
        );
        assert_eq!(
            born_cross_section_factors(&model, &photon_fusion, true).unwrap(),
            factors(&[-2, -2])
        );
    }

    #[test]
    fn coupling_budgets_preserve_beam_exchange_orbits() {
        let model = Model::qcd_qed();
        let generator = Generator::new(model.clone());
        let orders = |diagram: &FeynmanDiagram| {
            let mut total = BTreeMap::new();
            for (_, vertex) in diagram.vertices() {
                let rule = model
                    .vertex_rule_by_id(vertex.interaction.unwrap())
                    .unwrap();
                for (name, power) in rule.coupling_orders(&model) {
                    *total.entry(name).or_insert(0_usize) += power;
                }
            }
            total
        };
        let selected = |diagram: &FeynmanDiagram| orders(diagram).get("QCD") == Some(&2);
        // Filtering can renumber names; every other serialized field must agree,
        // including IDs, symmetry factors, projectors, fragments, cuts and routing.
        let snapshot = |diagram: &FeynmanDiagram| {
            diagram
                .clone()
                .with_name("budget-comparison")
                .to_json()
                .unwrap()
        };
        let mut nontrivial_orbit = false;
        for (process, loops) in [
            (Process::new(["u", "u"], ["u", "u"]), 1),
            (Process::new(["a", "a"], ["u", "u~"]), 2),
        ] {
            let mut exact_weights = None;
            for symmetrize_initial in [false, true] {
                let options = GenerationOptions::default()
                    .with_loop_count(loops, loops)
                    .unwrap()
                    .threads(1)
                    .symmetrize_initial(symmetrize_initial)
                    .symmetrize_left_right(false);
                let all = generator
                    .generate(&process, &options, GenerationType::CrossSection)
                    .unwrap();
                assert!(all.report.completed);
                let mut expected = all
                    .diagrams
                    .iter()
                    .filter(|diagram| selected(diagram))
                    .map(snapshot)
                    .collect::<Vec<_>>();
                assert!(!expected.is_empty(), "no QCD^2 diagrams at {loops} loops");
                assert!(
                    expected.len() < all.diagrams.len(),
                    "the budget must reject diagrams"
                );
                let filtered = generator
                    .generate(
                        &process,
                        &options.with_graph_filter(GenerationFilter::CouplingOrders(
                            BTreeMap::from([("QCD".to_owned(), (2, Some(2)))]),
                        )),
                        GenerationType::CrossSection,
                    )
                    .unwrap();
                assert!(filtered.report.completed);
                let mut actual = filtered.diagrams.iter().map(snapshot).collect::<Vec<_>>();
                expected.sort();
                actual.sort();
                assert_eq!(
                    actual, expected,
                    "{loops} loops, beam exchange {symmetrize_initial}"
                );

                // Compare each order class separately: signed totals across classes can cancel.
                let mut weights = BTreeMap::new();
                for diagram in &all.diagrams {
                    *weights.entry(orders(diagram)).or_insert_with(Atom::zero) +=
                        evaluated_overall_factor(diagram);
                }
                if symmetrize_initial {
                    assert_eq!(Some(weights), exact_weights, "{loops} loops");
                    // Two beams have only singleton and size-two orbits. Require the
                    // budgeted cases to exercise a size-two orbit, beyond fixed points.
                    nontrivial_orbit |= filtered.diagrams.iter().any(|diagram| {
                        let multiplicity = evaluated_overall_factor(diagram)
                            * Atom::num(i64::try_from(diagram.symmetry_factor()).unwrap());
                        multiplicity == Atom::num(2) || multiplicity == Atom::num(-2)
                    });
                } else {
                    exact_weights = Some(weights);
                }
            }
        }
        assert!(
            nontrivial_orbit,
            "the budget must retain a nontrivial beam orbit"
        );
    }

    #[test]
    fn cross_sections_reject_final_state_and_multi_beam_symmetrization() {
        let model = Model::qed();
        let options = GenerationOptions::default()
            .with_loop_count(1, 1)
            .unwrap()
            .threads(1);
        let generate = |process: Process, options: &GenerationOptions| {
            Generator::new(model.clone()).generate(&process, options, GenerationType::CrossSection)
        };
        assert!(matches!(
            generate(
                Process::new(["e-", "e+"], ["mu-", "mu+"]),
                &options.clone().symmetrize_final(true)
            ),
            Err(GenerationError::CrossSectionFinalStateSymmetrization)
        ));
        // Three identical beams, or two among three incoming particles, have no exchange that
        // a Lorentz transformation realizes.
        for incoming in [["a", "a", "a"], ["a", "a", "e-"]] {
            assert!(matches!(
                generate(
                    Process::new(incoming, ["e-", "e+"]),
                    &options.clone().symmetrize_initial(true)
                ),
                Err(GenerationError::UnsupportedBeamSymmetrization { incoming: 3 })
            ));
        }
    }

    #[test]
    fn left_right_mirror_images_merge_into_one_diagram() {
        let process = Process::new(["e-", "e+"], ["mu-", "mu+"]);
        let generate = |left_right| {
            Generator::new(Model::qed())
                .generate(
                    &process,
                    &GenerationOptions::default()
                        .with_loop_count(2, 2)
                        .unwrap()
                        .threads(1)
                        .symmetrize_left_right(left_right),
                    GenerationType::CrossSection,
                )
                .unwrap()
                .diagrams
        };
        let (exact, symmetrized) = (generate(false), generate(true));
        let ids = symmetrized
            .iter()
            .map(|diagram| diagram.id().unwrap())
            .collect::<BTreeSet<_>>();
        assert_eq!(ids.len(), symmetrized.len());
        assert!(symmetrized.len() < exact.len());
        assert_eq!(
            Atom::add_many(symmetrized.iter().map(evaluated_overall_factor)),
            Atom::add_many(exact.iter().map(evaluated_overall_factor))
        );
    }

    #[test]
    fn one_particle_irreducibility_keeps_one_point_loop_graphs() {
        // the only 1PI two-loop scalar tadpole: the leg vertex joined to both ends of a bubble
        let options = GenerationOptions::default()
            .with_loop_count(2, 2)
            .unwrap()
            .threads(1)
            .with_graph_filter(GenerationFilter::LoopOneParticleIrreducible);
        let result = Generator::new(scalar_model())
            .generate(
                &Process::new(["phi"], Vec::<&str>::new()),
                &options,
                GenerationType::Amplitude,
            )
            .unwrap();
        assert_eq!(result.diagrams.len(), 1);
        let diagram = &result.diagrams[0];
        assert_eq!(diagram.vertices().count(), 3);
        let mut lines = BTreeMap::<_, usize>::new();
        for (_, ends, edge) in diagram.edges() {
            if edge.external.is_none() {
                let (a, b) = (ends.source.unwrap(), ends.target.unwrap());
                *lines.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        let mut multiplicities = lines.into_values().collect::<Vec<_>>();
        multiplicities.sort();
        assert_eq!(multiplicities, [1, 1, 2]);
    }

    #[test]
    fn external_canonical_classes_keep_global_left_right_symmetry_out_of_raw_generation() {
        let model = fermion_model();
        let classes = |process: Process, generation_type, options: GenerationOptions| {
            ResolvedProcess::new(&model, &process, generation_type, &options)
                .unwrap()
                .externals(&model)
                .unwrap()
                .legs()
                .map(|(label, _)| label)
                .collect::<Vec<_>>()
        };
        let process = Process::new([1_i64, -1, 22], Vec::<i64>::new());
        assert_eq!(
            classes(
                process.clone(),
                GenerationType::Amplitude,
                GenerationOptions::default().symmetrize_initial(true)
            ),
            vec![
                ExternalLabel::Exact {
                    state: ExternalState::Incoming,
                    index: 0,
                },
                ExternalLabel::Exact {
                    state: ExternalState::Incoming,
                    index: 1,
                },
                ExternalLabel::WithinState(ExternalState::Incoming),
            ]
        );
        assert_eq!(
            classes(
                process,
                GenerationType::Amplitude,
                GenerationOptions::default()
                    .symmetrize_initial(true)
                    .symmetrize_external_fermions(true)
            ),
            vec![
                ExternalLabel::WithinState(ExternalState::Incoming),
                ExternalLabel::WithinState(ExternalState::Incoming),
                ExternalLabel::WithinState(ExternalState::Incoming),
            ]
        );

        // Sewn pairs are permuted after generation, so cross-section legs are always exact.
        let cross_section = Process::new([1_i64, -1], [22_i64]);
        assert_eq!(
            classes(
                cross_section,
                GenerationType::CrossSection,
                GenerationOptions::default().symmetrize_initial(true)
            ),
            vec![
                ExternalLabel::Exact {
                    state: ExternalState::Incoming,
                    index: 0,
                },
                ExternalLabel::Exact {
                    state: ExternalState::Incoming,
                    index: 1,
                },
                ExternalLabel::Exact {
                    state: ExternalState::Outgoing,
                    index: 2,
                },
                ExternalLabel::Exact {
                    state: ExternalState::Outgoing,
                    index: 3,
                },
            ]
        );

        let left_right = Process::new([1_i64, -1], [22_i64]);
        assert_eq!(
            classes(
                left_right,
                GenerationType::CrossSection,
                GenerationOptions::default().symmetrize_left_right(true)
            ),
            vec![
                ExternalLabel::Exact {
                    state: ExternalState::Incoming,
                    index: 0,
                },
                ExternalLabel::Exact {
                    state: ExternalState::Incoming,
                    index: 1,
                },
                ExternalLabel::Exact {
                    state: ExternalState::Outgoing,
                    index: 2,
                },
                ExternalLabel::Exact {
                    state: ExternalState::Outgoing,
                    index: 3,
                },
            ]
        );
    }

    #[test]
    fn numerator_construction_localizes_internal_ufo_data() {
        let source = fermion_model().to_json().unwrap();
        let mut definition = serde_json::from_str::<serde_json::Value>(&source).unwrap();
        let fermion = definition["particles"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|particle| particle["name"] == "f")
            .unwrap();
        fermion["propagator"] = serde_json::json!("custom_f_prop");
        definition["propagators"]
            .as_array_mut()
            .unwrap()
            .push(serde_json::json!({
                "name": "custom_f_prop",
                "particle": "f",
                "numerator": "UFO::{}::P(UFO::{}::idx(1,1))*UFO::{}::P(UFO::{}::idx(1,2))*UFO::{}::Metric(UFO::{}::idx(1,1),UFO::{}::dummy(1))*UFO::{}::Metric(UFO::{}::dummy(1),UFO::{}::idx(1,2))",
                "denominator": "P^2"
            }));
        let definition = serde_json::to_string(&definition).unwrap();
        let generator = Generator::new(Model::from_json(&definition).unwrap());

        let mut graph = Graph::new();
        let interaction = graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V")));
        let external = graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Outgoing,
            particle: particle(&generator.model, 22),
        }));
        graph
            .add_edge(
                interaction,
                interaction,
                true,
                particle(&generator.model, 1),
            )
            .unwrap();
        graph
            .add_edge(interaction, external, false, particle(&generator.model, 22))
            .unwrap();

        let diagram = generator
            .to_diagram(
                Some("localized".to_owned()),
                GenerationType::Amplitude,
                ColoredTopology {
                    graph: graph.build(),
                    symmetry: 1,
                    cut_partitions: Vec::new(),
                    topology_threshold_candidates: Vec::new(),
                },
                FermionSigns {
                    internal_loop_negative: true,
                    ..FermionSigns::default()
                },
                1,
                &GenerationOptions::default(),
            )
            .unwrap();
        let internal_numerator = diagram
            .edges()
            .find(|(_, _, edge)| edge.external.is_none())
            .unwrap()
            .2
            .numerator
            .to_plain_string();
        assert!(internal_numerator.contains("gammalooprs::Q"));
        assert!(internal_numerator.contains("gammalooprs::hedge"));
        assert!(internal_numerator.contains("gammalooprs::hedge"));
        assert!(internal_numerator.contains("gammalooprs::edge"));
        assert!(!internal_numerator.contains("idx("));
        assert!(!internal_numerator.contains("Slash(P)"));

        let momentum = symbol!("gammalooprs::Q");
        assert!(momentum.has_tag("spenso::tensor"));
        assert!(momentum.has_tag("spenso::rank1"));

        let external_numerator = &diagram
            .edges()
            .find(|(_, _, edge)| edge.external.is_some())
            .unwrap()
            .2
            .numerator;
        assert_eq!(external_numerator, &Atom::one());

        let vertex_numerator = diagram
            .vertex(feynkit_graph::VertexId(interaction))
            .unwrap()
            .numerator
            .to_plain_string();
        assert!(vertex_numerator.contains("gammalooprs::hedge"));
        assert!(vertex_numerator.contains("gammalooprs::hedge"));
    }
}
