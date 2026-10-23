use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt;
use std::hash::{Hash, Hasher};
use std::sync::{Arc, Mutex, RwLock};

use feynkit_graph::symbols;
use feynkit_graph::{
    DiagramCut, DiagramCutSide, DiagramEdge, DiagramEndpoint, DiagramError, DiagramHalfEdge,
    DiagramId, DiagramThresholdCandidate, DiagramVertex, EdgeId, ExternalState, FeynmanDiagram,
    NumeratorInterface, VertexSlot,
};
use feynkit_model::{Model, ModelError, Particle, ParticleId, VertexRule, VertexRuleId};
use idenso::color::CS;
use idenso::dirac::AGS;
use idenso::epsilon::EPSILON_SYMBOL;
use idenso::representations::Bispinor;
use idenso::tensor::SymbolicTensor;
use idenso::tensor::composition;
use idenso::tensor::inference::InterfaceInference;
use idenso::{CookMode, CookSettings};
use linnet::half_edge::involution::Flow;
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike};
use linnet::half_edge::{HedgeGraph, NodeIndex};
use linnet::permutation::Permutation;
use rayon::ThreadPoolBuildError;
#[cfg(not(target_arch = "wasm32"))]
use rayon::ThreadPoolBuilder;
#[cfg(not(target_arch = "wasm32"))]
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};
use serde::{Deserialize, Serialize};
use spenso::network::library::symbolic::ETS;
use spenso::structure::abstract_index::{AIND_SYMBOLS, AbstractIndex};
use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
use spenso::structure::slot::{IsAbstractSlot, ParseableAind, SlotMatch, SlotMatcher};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, UserData};
use symbolica::graph::{GenerationSettings, Graph, HalfEdge};
use symbolica::symbol;
use thiserror::Error;

use feynkit_amplitude::ColorRepresentation;

use crate::{
    FilterScope, GenerationFilter, GenerationFilterKind, GenerationOptions, GenerationType,
    GroupingError, NumeratorGrouping, ParticleSelector, Process, ProcessError, SelectorError,
    SelfEnergyFilterOptions, SewnFilterOptions, SnailFilterOptions, TadpoleFilterOptions,
    VertexSelector, grouping, momentum_symbol,
};

#[cfg(test)]
#[path = "main_generation_regressions.rs"]
mod main_generation_regressions;

#[derive(Debug, Error)]
pub enum GenerationError {
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
    #[error("graph automorphism factor cannot be represented by u64: {0}")]
    SymmetryFactor(String),
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
    #[error("numerator grouping requires canonical comparison or at least one numerical sample")]
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
    #[error("spanning-tree edge between nodes {node} and {parent} is missing")]
    MissingSpanningTreeEdge { node: usize, parent: usize },
    #[error("external node {node} has no incident edge")]
    MissingExternalEdge { node: usize },
    #[error("external leg {index} has no resolved particle")]
    MissingExternalParticle { index: usize },
    #[error("projector refers to unknown process external leg {index}")]
    UnknownProjectorExternalIndex { index: usize },
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

#[derive(Debug, Clone)]
pub struct GroupMember {
    /// Stable generated-order index assigned before zero-numerator removal.
    pub source_diagram: usize,
    /// Stable content-derived ID of the source diagram.
    pub source_id: DiagramId,
    /// Finalized display name of the source diagram.
    pub source_name: String,
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
    pub topology_count: usize,
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

struct NumeratorInstantiation<'a> {
    owner: NumeratorOwner,
    legs: &'a [NumeratorHalfEdge],
}

impl NumeratorInstantiation<'_> {
    fn propagator_expression(
        &self,
        propagator: &Atom,
        color: i64,
    ) -> Result<Atom, GenerationError> {
        let spin = self.instantiate(propagator, NumeratorSector::Spin)?;
        Ok(match self.edge_color_identity(color)? {
            None => spin,
            Some(color) => spin * color,
        })
    }

    /// Bind canonical template labels simultaneously. Momentum IDs may alias;
    /// half-edge and dummy-owner labels retain their distinct identities.
    fn bind_template(&self, expression: &Atom) -> Result<Atom, GenerationError> {
        let mut error = None;
        let expression = expression.replace_map_bottom_up(|term, _, out| {
            let AtomView::Fun(function) = term else {
                return;
            };
            let head = function.get_symbol();
            let mut arguments = function.iter();
            let owner = if head == symbols::hedge_index() || head == symbols::momentum() {
                let leg = usize::try_from(
                    arguments
                        .next()
                        .expect("localized template momentum and half-edge heads carry a leg ID"),
                )
                .expect("localized template leg IDs are nonnegative integers");
                if head == symbols::momentum() {
                    self.legs[leg].edge
                } else {
                    self.legs[leg].hedge
                }
            } else if head == symbols::vertex_index() || head == symbols::edge_index() {
                arguments.next();
                match self.owner {
                    NumeratorOwner::Vertex(owner) | NumeratorOwner::Edge(owner) => owner,
                }
            } else {
                return;
            };
            match i64::try_from(owner) {
                Ok(owner) => {
                    **out = FunctionBuilder::new(head)
                        .add_arg(owner)
                        .add_args(arguments)
                        .finish();
                }
                Err(_) => {
                    error = Some(GenerationError::ArithmeticOverflow(
                        "a symbolic numerator owner",
                    ));
                }
            }
        });
        error.map_or(Ok(expression), Err)
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
                // UFO's real Levi-Civita tensor has ε^{0123} = −1 (ALOHA,
                // matching FeynRules), so Tr[γ5 γ^μ γ^ν γ^ρ γ^σ] = 4i ε^{μνρσ}
                // for the same γ5 = diag(−1,−1,1,1). Idenso's epsilon is
                // normalized by Tr[…] = 4 ε^{μνρσ}, i.e. ε^{0123} = −i.
                "Epsilon" => self.exact_arguments(term, &arguments, 4).map(|args| {
                    let epsilon = args
                        .iter()
                        .fold(
                            FunctionBuilder::new(*EPSILON_SYMBOL),
                            |builder, argument| builder.add_arg(self.minkowski(argument.clone())),
                        )
                        .finish();
                    -Atom::i() * epsilon
                }),
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
            if function.get_symbol() != momentum_symbol() || function.get_nargs() != 2 {
                return;
            }
            let mut arguments = function.iter();
            let edge = arguments.next().expect("the arity was checked");
            let index = arguments.next().expect("the arity was checked");
            if Self::is_representation(index, "mink") {
                return;
            }
            **out = FunctionBuilder::new(momentum_symbol())
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
        let dummy_symbol = match self.owner {
            NumeratorOwner::Vertex(_) => symbols::vertex_index(),
            NumeratorOwner::Edge(_) => symbols::edge_index(),
        };
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
            let expected = match self.owner {
                NumeratorOwner::Vertex(_) => symbols::vertex_index(),
                NumeratorOwner::Edge(_) => symbols::edge_index(),
            };
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
        if function.get_symbol() != momentum_symbol() || function.get_nargs() != 1 {
            return Err(self.invalid_tensor(momentum, "expected a rank-one localized momentum"));
        }
        Ok(FunctionBuilder::new(momentum_symbol())
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
            if [symbols::vertex_index(), symbols::edge_index()]
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
        Ok(if shift == 0 {
            symbols::hedge_index().call((hedge,))
        } else {
            symbols::hedge_index().call((hedge, shift))
        })
    }

    fn dummy(&self, local: usize) -> Result<Atom, GenerationError> {
        let (symbol, owner) = match self.owner {
            NumeratorOwner::Vertex(vertex) => (symbols::vertex_index(), vertex),
            NumeratorOwner::Edge(edge) => (symbols::edge_index(), edge),
        };
        let owner = i64::try_from(owner)
            .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic numerator owner"))?;
        let local = i64::try_from(local)
            .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic dummy index"))?;
        Ok(if local == 0 {
            symbol.call((owner,))
        } else {
            symbol.call((owner, local))
        })
    }

    fn momentum(
        &self,
        leg: NumeratorHalfEdge,
        index: Option<Atom>,
    ) -> Result<Atom, GenerationError> {
        let edge = i64::try_from(leg.edge)
            .map_err(|_| GenerationError::ArithmeticOverflow("a symbolic edge identifier"))?;
        let mut momentum = FunctionBuilder::new(momentum_symbol()).add_arg(edge);
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
pub(crate) struct Generator {
    model: Arc<Model>,
    numerators: Arc<NumeratorCache>,
}

/// Localized factors depend on the immutable model and their index assignment,
/// not on the rest of the diagram. Keep them for this generator's lifetime.
#[derive(Default)]
struct NumeratorCache {
    interactions: RwLock<HashMap<InteractionNumeratorKey, SymbolicTensor<PartialStructure>>>,
    propagators: RwLock<HashMap<PropagatorNumeratorKey, SymbolicTensor<PartialStructure>>>,
    interaction_templates: RwLock<HashMap<InteractionNumeratorKey, Arc<NumeratorTemplate>>>,
    propagator_templates: RwLock<HashMap<ParticleId, Arc<NumeratorTemplate>>>,
}

/// The same model-bound owner retains a rule's admission independently of the
/// graph's owner IDs. Dynamic source notation must be lowered at its actual owner.
enum NumeratorTemplate {
    Dynamic,
    Admitted {
        tensor: Box<SymbolicTensor<PartialStructure>>,
        labels: Vec<(AbstractIndex, Atom)>,
    },
}

impl NumeratorTemplate {
    fn supports(expression: &Atom) -> bool {
        if !InterfaceInference::normalization_is_intrinsic(expression.as_view()) {
            return false;
        }
        let mut supported = true;
        let mut slots = SlotMatcher::default();
        expression.visitor(&mut |node| {
            // Admission may materialize these occurrences into fresh dummies,
            // even when a surrounding power or chain has no public ports.
            // Keep their allocation in the actual owner's construction.
            supported &= slots.compact_representation(node).is_none();
            if let AtomView::Fun(function) = node {
                // The occurrence-order owner cannot enumerate ports hidden in
                // bundles. Reject them before an additive term can move first
                // under owner-dependent normalization.
                supported &= function.get_symbol() != AIND_SYMBOLS.aind;
                // These labels must originate in localization. Literal graph
                // IDs could capture a generated label in the canonical frame.
                supported &= ![
                    symbols::hedge_index(),
                    symbols::vertex_index(),
                    symbols::edge_index(),
                    symbols::momentum(),
                    AIND_SYMBOLS.openind,
                ]
                .contains(&function.get_symbol());
            }
            if let SlotMatch::Explicit(slot) = slots.classify(node) {
                supported &= InterfaceInference::normalization_is_intrinsic(slot.index());
            }
            if let AtomView::Var(variable) = node
                && let UserData::Atom(payload) = variable.get_symbol().get_data()
            {
                supported &= Self::supports(payload);
            }
            supported
        });
        supported
    }

    fn admitted(tensor: SymbolicTensor<PartialStructure>) -> Result<Self, GenerationError> {
        let tensor = tensor.finish_parts().map_err(DiagramError::from)?;
        // Anonymous occurrences need the checked positional correspondence
        // owner. Intrinsic explicit templates have a complete label bijection.
        if !tensor.structure().open_positions().is_empty()
            || !InterfaceInference::normalization_is_intrinsic(tensor.expression().as_view())
            || composition::interface_in_atom_order(tensor.expression(), tensor.structure())
                .is_err()
        {
            return Ok(Self::Dynamic);
        }
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let mut indices = HashSet::new();
        let mut slots = SlotMatcher::default();
        let mut fixed_dimensions = true;
        let mut collectable_indices = true;
        tensor.expression().visitor(&mut |node| {
            if let SlotMatch::Explicit(slot) = slots.classify(node) {
                if !slot.is_concrete_index() {
                    match AbstractIndex::try_from(slot.index()) {
                        Ok(index) => {
                            indices.insert(index);
                        }
                        // Slot-shaped opaque scalar metadata need not contain
                        // an abstract index, even in an admitted expression.
                        Err(_) => collectable_indices = false,
                    }
                }
                fixed_dimensions &= Self::supports(&cooking.uncook(slot.dimension()));
            }
            collectable_indices
        });
        if !collectable_indices {
            return Ok(Self::Dynamic);
        }
        for slot in tensor.structure().logical_slots() {
            if let PartialIndex::Explicit(index) = slot.aind {
                indices.insert(index);
            }
        }
        let labels = indices
            .into_iter()
            .map(|index| (index, cooking.uncook(index.to_atom().as_view())))
            .collect::<Vec<_>>();
        let mut fixed_labels = true;
        for (_, label) in &labels {
            fixed_labels &= InterfaceInference::normalization_is_intrinsic(label.as_view());
            label.visitor(&mut |node| {
                // An opaque index depending on a momentum ID could itself
                // alias when incident legs share an edge.
                fixed_labels &= node.get_symbol() != Some(symbols::momentum());
                fixed_labels
            });
        }
        if !fixed_dimensions || !fixed_labels {
            return Ok(Self::Dynamic);
        }
        let public_keys = tensor
            .structure()
            .logical_slots()
            .into_iter()
            .filter_map(|slot| match slot.aind {
                PartialIndex::Explicit(index) => {
                    Some(slot.rep().base().slot::<AbstractIndex, _>(index))
                }
                PartialIndex::Open(_) => None,
            })
            .collect::<HashSet<_>>();
        let representations = public_keys
            .iter()
            .map(|slot| slot.rep().rep.base())
            .collect::<HashSet<_>>();
        let mut scoped_alias = false;
        tensor.expression().visitor(&mut |scope| {
            if matches!(scope, AtomView::Add(_) | AtomView::Pow(_)) {
                for &representation in &representations {
                    let paired = InterfaceInference::has_explicit_index_pairs(
                        scope,
                        representation,
                        |leaf| {
                            let AtomView::Fun(function) = leaf else {
                                return true;
                            };
                            !function.iter().any(|argument| {
                                slots
                                    .parse::<LibraryRep, AbstractIndex>(argument)
                                    .is_ok_and(|slot| {
                                        let key =
                                            slot.rep().base().slot::<AbstractIndex, _>(slot.aind());
                                        public_keys.contains(&key)
                                    })
                            })
                        },
                    );
                    // A local bound pair must not be mistaken for a later
                    // public occurrence by the source-order helper. Multi-port
                    // leaves and detached chain scopes can be inconclusive;
                    // retain actual-owner admission in those cases too.
                    if !matches!(paired, Ok(false)) {
                        scoped_alias = true;
                        return false;
                    }
                }
            }
            !scoped_alias
        });
        if scoped_alias {
            return Ok(Self::Dynamic);
        }
        // Complete the source metadata before the retained template is cloned.
        // The identity map keeps its established interface and collision history.
        let tensor = tensor
            .with_renamed_indices(tensor.expression().clone(), &HashMap::new())
            .map_err(DiagramError::from)?;
        Ok(Self::Admitted {
            tensor: Box::new(tensor),
            labels,
        })
    }

    fn instantiate(
        &self,
        instantiation: &NumeratorInstantiation<'_>,
    ) -> Result<SymbolicTensor<PartialStructure>, GenerationError> {
        let Self::Admitted { tensor, labels } = self else {
            unreachable!("dynamic rules are constructed at their actual owner")
        };
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let mut renaming = HashMap::with_capacity(labels.len());
        let mut encoded = HashMap::new();
        for (index, decoded) in labels {
            let mut mapped = instantiation.bind_template(decoded)?;
            let target = if mapped == decoded {
                *index
            } else {
                let mut scopes = Vec::new();
                while let Some(function) = mapped.as_fun_view()
                    && function.get_symbol() == AbstractIndex::scope_symbol()
                {
                    let arguments = function.iter().collect::<Vec<_>>();
                    let [AtomView::Var(scope), inner] = arguments.as_slice() else {
                        unreachable!("admitted index scope")
                    };
                    scopes.push(scope.get_symbol());
                    mapped = (*inner).to_owned();
                }
                if let Some(function) = mapped.as_fun_view() {
                    mapped = cooking.cook_function(function).map_err(|source| {
                        DiagramError::TensorCooking {
                            operation: "instantiating an admitted numerator template",
                            source: Box::new(source),
                        }
                    })?;
                }
                scopes.into_iter().rev().fold(
                    AbstractIndex::try_from(mapped.as_view()).map_err(|source| {
                        DiagramError::TensorIndex {
                            operation: "instantiating an admitted numerator template",
                            source: Box::new(source),
                        }
                    })?,
                    |index, scope| index.scoped(scope),
                )
            };
            renaming.insert(*index, target);
            // Reversible labels are opaque during the native-label pass. Map
            // each complete index afterwards, including symbols inside scopes.
            if decoded != &index.to_atom() {
                encoded.insert(index.to_atom(), target.to_atom());
            }
        }
        let expression = instantiation.bind_template(tensor.expression())?;
        let expression = if encoded.is_empty() {
            expression
        } else {
            expression.replace_map(|node, _, output| {
                if let Some(target) = encoded.get(&node.to_owned()) {
                    **output = target.clone();
                }
            })
        };
        // Admission established the raw rule; scalar momentum aliases retain
        // its metadata. This injective relabeling also transports identities
        // retained only by a zero. No callbacks run.
        Ok(tensor
            .with_renamed_indices(expression, &renaming)
            .map_err(DiagramError::from)?)
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
    /// The diagram's routing, overall factor, numerator prefactor and projector
    /// are retained. Existing local numerator fragments are replaced with the
    /// model's Feynman rules, just as during diagram generation.
    fn apply_feynman_rules(&self) -> Result<FeynmanDiagram, GenerationError>;

    /// Construct the external wavefunctions in the finalized half-edge frame.
    /// Callers with an explicit projector retain that supplied expression.
    fn external_state_projector(&self)
    -> Result<SymbolicTensor<PartialStructure>, GenerationError>;
}

impl FeynmanDiagramRulesExt for FeynmanDiagram {
    fn apply_feynman_rules(&self) -> Result<FeynmanDiagram, GenerationError> {
        self.validate()?;
        Generator::new(self.model_arc()).construct_numerators(self.clone())
    }

    fn external_state_projector(
        &self,
    ) -> Result<SymbolicTensor<PartialStructure>, GenerationError> {
        let generator = Generator::new(self.model_arc());
        let mut projector = Atom::one();
        for (edge, endpoints, data) in self.edges() {
            let Some(external) = &data.external else {
                continue;
            };
            for (vertex, endpoint, state) in [
                (
                    endpoints.source,
                    DiagramEndpoint::Source,
                    ExternalState::Outgoing,
                ),
                (
                    endpoints.target,
                    DiagramEndpoint::Target,
                    ExternalState::Incoming,
                ),
            ] {
                if vertex.is_none() {
                    continue;
                }
                let state = if endpoints.source.is_some() && endpoints.target.is_some() {
                    state
                } else {
                    external.state
                };
                let hedge = self
                    .half_edge_id(DiagramHalfEdge { edge, endpoint })
                    .expect("finalized external endpoint")
                    .0;
                projector *=
                    generator.external_wavefunction(state, data.particle, edge.0, hedge)?;
            }
        }
        Ok(NumeratorInterface::admit(projector, None)?)
    }
}

impl Generator {
    pub fn new(model: impl Into<Arc<Model>>) -> Self {
        Self {
            model: model.into(),
            numerators: Arc::default(),
        }
    }

    pub fn generate(
        &self,
        process: &Process,
        options: &GenerationOptions,
        generation_type: GenerationType,
    ) -> Result<GenerationResult, GenerationError> {
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
        let options = &options;
        self.validate_options(generation_type, options)?;
        process.validate_covariant_cut_filters(&self.model, options, generation_type)?;
        let expanded_process = process.clone().with_final_state_alternatives(
            process.covariant_cut_states(&self.model, generation_type)?,
        )?;
        let process = &expanded_process;
        let resolved = ResolvedProcess::new(&self.model, process, generation_type, options)?;
        let signatures = InteractionSignatures::new(&self.model, options)?;
        let external_edges = resolved.external_edges(&self.model)?;
        let generation_signatures: Vec<_> = signatures.keys().cloned().collect();

        let progress_count = Arc::new(Mutex::new(0));
        let progress_options = options.clone();
        options.report_progress("topologies", 0, None);
        // Fast cut filter switch multiplicity is no longer a separate CLI option;
        // FeynKit owns the physical-cut filtering strategy.
        let mut settings = GenerationSettings::new()
            .max_loops(*options.loop_count().end())
            .allow_self_loops(options.allow_self_loops)
            .allow_zero_flow_edges(options.allow_zero_flow_edges)
            .progress_fn(Box::new(move |_| {
                let mut count = progress_count.lock().unwrap();
                *count += 1;
                progress_options.report_progress("topologies", *count, None);
                !progress_options.cancellation_requested()
            }));
        if let Some(filter) = options.filter.clone() {
            let filter_options = options.clone();
            settings = settings.filter_fn(Box::new(move |graph, completed_vertices| {
                !filter_options.cancellation_requested() && filter(graph, completed_vertices)
            }));
        }
        if let Some(maximum) = options.max_vertices {
            // Symbolica counts the degree-one external nodes as vertices. The
            // public FeynKit option counts interaction vertices instead.
            settings = settings.max_vertices(maximum.checked_add(external_edges.len()).ok_or(
                GenerationError::ArithmeticOverflow("the Symbolica vertex limit"),
            )?);
        }
        if let Some(maximum) = options.max_bridges(*options.loop_count().start()) {
            settings = settings.max_bridges(maximum);
        }
        let abort_options = options.clone();
        settings = settings.abort_check(Box::new(move || abort_options.cancellation_requested()));

        let generated = Graph::generate(&external_edges, &generation_signatures, settings);
        let (raw_graphs, mut completed) = match generated {
            Ok(graphs) => (graphs, true),
            Err(graphs) => (graphs, false),
        };
        if options.cancellation_requested() {
            completed = false;
        }

        let mut raw_graphs: Vec<_> = raw_graphs.into_iter().collect();
        // Sort for reproducible generated IDs, then discard contributions below
        // the requested loop range before more expensive processing.
        raw_graphs.sort_by(|left, right| left.0.cmp(&right.0));
        raw_graphs.retain(|(graph, _)| options.loop_count().contains(&graph.num_loops()));
        let mut filtered_graphs = Vec::with_capacity(raw_graphs.len());
        let total = raw_graphs.len();
        options.report_progress("topology_filters", 0, Some(total));
        for (index, entry) in raw_graphs.into_iter().enumerate() {
            if options.cancellation_requested() {
                completed = false;
                break;
            }
            if self.passes_topology_filters(&entry.0, options)? {
                filtered_graphs.push(entry);
            }
            options.report_progress("topology_filters", index + 1, Some(total));
        }
        let raw_graphs = filtered_graphs;
        let topology_count = raw_graphs.len();

        let assignments = Mutex::new(0);
        options.report_progress("interactions", 0, Some(topology_count));
        let color_graph = |(graph, symmetry): &(
            Graph<NodeColor, EdgeColor>,
            symbolica::domains::integer::Integer,
        )| {
            if options.cancellation_requested() {
                return Ok(None);
            }
            let symmetry = match u64::try_from(symmetry.clone()) {
                Ok(symmetry) => symmetry,
                Err(error) => {
                    return Err(GenerationError::SymmetryFactor(error.to_string()));
                }
            };
            let result = self
                .assign_interactions(graph, &signatures)
                .map(|graphs| (graphs, symmetry));
            let mut count = assignments.lock().unwrap();
            *count += 1;
            options.report_progress("interactions", *count, Some(topology_count));
            result.map(Some)
        };
        // An explicit single thread runs every stage on the caller without a
        // Rayon pool. Otherwise respect an explicit thread count, or Rayon's
        // default logical-CPU count. Browser kernels process the same graphs serially.
        #[cfg(not(target_arch = "wasm32"))]
        let pool = if options.threads == Some(1) {
            None
        } else {
            let mut pool = ThreadPoolBuilder::new();
            if let Some(threads) = options.threads {
                pool = pool.num_threads(threads);
            }
            Some(pool.build()?)
        };
        #[cfg(target_arch = "wasm32")]
        let colored = raw_graphs
            .iter()
            .map(color_graph)
            .collect::<Result<Vec<_>, _>>()?;
        #[cfg(not(target_arch = "wasm32"))]
        let colored = if let Some(pool) = &pool {
            pool.install(|| {
                raw_graphs
                    .par_iter()
                    .map(color_graph)
                    .collect::<Result<Vec<_>, _>>()
            })?
        } else {
            raw_graphs
                .iter()
                .map(color_graph)
                .collect::<Result<Vec<_>, _>>()?
        };
        if options.cancellation_requested() {
            completed = false;
        }
        let mut colored: Vec<_> = colored
            .into_iter()
            .flatten()
            .flat_map(|(graphs, symmetry)| {
                graphs
                    .into_iter()
                    .map(move |(graph, multiplicity)| ColoredTopology {
                        graph,
                        symmetry,
                        multiplicity,
                        cut_partitions: Vec::new(),
                        topology_threshold_candidates: Vec::new(),
                    })
            })
            .collect();
        colored.sort_by(|left, right| left.graph.cmp(&right.graph));
        let interaction_assignment_count = colored.len();

        let mut diagram_inputs = Vec::new();
        options.report_progress("interaction_filters", 0, Some(interaction_assignment_count));
        for (index, mut topology) in colored.into_iter().enumerate() {
            options.report_progress(
                "interaction_filters",
                index,
                Some(interaction_assignment_count),
            );
            if options.cancellation_requested() {
                completed = false;
                break;
            }
            topology.graph = resolved.canonicalize_numerator_graph(
                &topology.graph,
                &self.model,
                &options.numerator_grouping,
            )?;
            if generation_type == GenerationType::CrossSection {
                let inventory = resolved.cut_partitions(&topology.graph, &self.model, options)?;
                topology.cut_partitions = inventory.physical;
                topology.topology_threshold_candidates = inventory.topology_threshold_candidates;
                if topology.cut_partitions.is_empty()
                    || !self.retain_cut_amplitude_partitions(&mut topology, options)?
                    || !resolved.passes_sewn_filter(&topology.graph, options)
                {
                    continue;
                }
            }
            topology.graph =
                ResolvedProcess::repair_cp_interaction_rules(topology.graph, &self.model)?;
            let fermion_loop_count = self.closed_fermion_loop_count(&topology.graph, false)?;
            if !self.passes_interaction_filters(&topology.graph, fermion_loop_count, options)? {
                continue;
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
                symmetry: topology.symmetry,
                multiplicity: topology.multiplicity,
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
            // already the representative. Its orientation-independent key also
            // provides the same deterministic ordering.
            let (topology, reversed_edges) =
                if matches!(options.numerator_grouping, NumeratorGrouping::None) {
                    (representative.clone(), BTreeSet::new())
                } else {
                    (topology, normalized.reversed_edges)
                };
            diagram_inputs.push((
                topology,
                representative,
                grassmann_loop_count,
                normalized.signs,
                reversed_edges,
            ));
        }
        if !options.cancellation_requested() {
            options.report_progress(
                "interaction_filters",
                interaction_assignment_count,
                Some(interaction_assignment_count),
            );
        }
        let total = diagram_inputs.len();
        options.report_progress("numerators", 0, Some(total));
        let width = diagram_inputs.len().saturating_sub(1).to_string().len();
        let constructed = Mutex::new(0);
        let construct = |(
            index,
            (comparison, representative, grassmann_loop_count, signs, reversed_edges),
        ): (
            usize,
            (
                ColoredTopology,
                ColoredTopology,
                usize,
                FermionSigns,
                BTreeSet<EdgeId>,
            ),
        )| {
            if options.cancellation_requested() {
                return Ok::<_, GenerationError>(None);
            }
            let name = format!("{}{index:0width$}", options.graph_prefix);
            let identical = comparison == representative;
            let signs = FermionSigns {
                internal_loop_negative: grassmann_loop_count % 2 == 1,
                ..signs
            };
            let internal_reversals = comparison
                .graph
                .edges()
                .iter()
                .enumerate()
                .filter(|(_, edge)| {
                    [edge.vertices.0, edge.vertices.1].iter().all(|vertex| {
                        matches!(
                            comparison.graph.nodes()[*vertex].data,
                            ColoredNode::Interaction(_)
                        )
                    })
                })
                .map(|(edge, _)| reversed_edges.contains(&EdgeId(edge)))
                .collect::<Vec<_>>();
            let comparison =
                self.to_diagram(name.clone(), generation_type, comparison, signs, options)?;
            let representative = if identical {
                comparison.clone()
            } else {
                self.to_diagram(name, generation_type, representative, signs, options)?
            };
            // Finalization orders external carriers before the internal edges.
            // Transport the comparison-to-representative reversals into that
            // finalized edge frame; process flow fixes all external carriers.
            let reversed_edges = comparison
                .edges()
                .filter(|(_, _, edge)| edge.external.is_none())
                .zip(internal_reversals)
                .filter_map(|((edge, _, _), reversed)| reversed.then_some(edge))
                .collect();
            let pair = DiagramPair {
                comparison,
                representative,
                reversed_edges,
            };
            let mut count = constructed.lock().unwrap();
            *count += 1;
            options.report_progress("numerators", *count, Some(total));
            Ok(Some(pair))
        };
        #[cfg(target_arch = "wasm32")]
        let diagram_pairs = diagram_inputs
            .into_iter()
            .enumerate()
            .map(construct)
            .collect::<Result<Vec<_>, _>>()?;
        #[cfg(not(target_arch = "wasm32"))]
        let diagram_pairs = if let Some(pool) = &pool {
            pool.install(|| {
                diagram_inputs
                    .into_par_iter()
                    .enumerate()
                    .map(construct)
                    .collect::<Result<Vec<_>, _>>()
            })?
        } else {
            diagram_inputs
                .into_iter()
                .enumerate()
                .map(construct)
                .collect::<Result<Vec<_>, _>>()?
        };
        let mut diagram_pairs = diagram_pairs.into_iter().flatten().collect::<Vec<_>>();
        if options.cancellation_requested() {
            completed = false;
        }
        options.report_progress("selection", 0, None);
        if !options.forced_cuts.is_empty() {
            let mut matched = vec![false; options.forced_cuts.len()];
            let mut retained = Vec::new();
            for mut pair in diagram_pairs {
                let cuts = pair
                    .comparison
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
                    pair.comparison = pair.comparison.with_cuts(cuts)?;
                    let partitions =
                        transport_cut_metadata(pair.comparison.cuts(), &pair.reversed_edges)
                            .into_iter()
                            .map(|cut| (cut.left.half_edges, cut.right.half_edges))
                            .collect();
                    pair.representative = pair.representative.with_cut_partitions(partitions)?;
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
        let mut keyed_diagrams = diagram_pairs
            .into_iter()
            .map(|pair| Ok((pair.comparison.canonical_key()?, pair)))
            .collect::<Result<Vec<_>, GenerationError>>()?;
        keyed_diagrams.sort_by(|left, right| left.0.cmp(&right.0));
        let width = keyed_diagrams
            .len()
            .saturating_sub(1)
            .to_string()
            .len()
            .max(1);
        let mut comparison_diagrams = Vec::with_capacity(keyed_diagrams.len());
        let mut representatives = Vec::with_capacity(keyed_diagrams.len());
        for (index, (_, mut pair)) in keyed_diagrams.into_iter().enumerate() {
            let name = format!("{}{index:0width$}", options.graph_prefix);
            let selected_by_id = options
                .selected_diagram_ids
                .as_ref()
                .is_some_and(|selected| {
                    selected.contains(&pair.comparison.id())
                        || selected.contains(&pair.representative.id())
                });
            let selected_by_name = options
                .selected_diagram_names
                .as_ref()
                .is_some_and(|selected| selected.contains(&name));
            let has_selection =
                options.selected_diagram_ids.is_some() || options.selected_diagram_names.is_some();
            if (has_selection && !selected_by_id && !selected_by_name)
                || options.vetoed_diagram_ids.contains(&pair.comparison.id())
                || options
                    .vetoed_diagram_ids
                    .contains(&pair.representative.id())
                || options.vetoed_diagram_names.contains(&name)
            {
                continue;
            }
            pair.comparison = pair.comparison.with_name(&name);
            pair.representative = pair.representative.with_name(&name);
            let requested_basis = options
                .loop_momentum_bases
                .get(&pair.representative.id())
                .or_else(|| options.loop_momentum_bases.get(&pair.comparison.id()))
                .or_else(|| options.named_loop_momentum_bases.get(&name));
            if let Some(edges) = requested_basis {
                pair.comparison = pair.comparison.with_loop_momentum_edges(edges)?;
                pair.representative = pair.representative.with_loop_momentum_edges(edges)?;
            }
            comparison_diagrams.push(pair.comparison);
            representatives.push((pair.representative, pair.reversed_edges));
        }
        let grouping = if options.cancellation_requested() {
            completed = false;
            &NumeratorGrouping::None
        } else {
            &options.numerator_grouping
        };
        let group = || {
            grouping::group_diagrams(
                comparison_diagrams,
                &self.model,
                grouping,
                options.symmetrizes_left_right(),
                options,
            )
        };
        #[cfg(target_arch = "wasm32")]
        let grouped = group()?;
        #[cfg(not(target_arch = "wasm32"))]
        let grouped = if let Some(pool) = &pool {
            pool.install(group)?
        } else {
            group()?
        };
        if options.cancellation_requested() {
            completed = false;
        }

        // Numerator-aware grouping must use the canonical comparison graphs:
        // normalizing their charged flows changes the momentum convention and
        // prevents mirror representatives from collapsing.  Runtime exposure,
        // however, uses the normalized representative selected by each group.
        // Transport the merged canonical cut inventories into that
        // representative's orientation and retain the comparison-derived
        // grouping factor and routing choice.
        let zero_numerator_count = grouped.zero_numerator_count;
        let mut groups = grouped.groups;
        let total = groups.len();
        options.report_progress_with_zero_count("grouping", 0, Some(total), zero_numerator_count);
        let mut diagrams = Vec::with_capacity(grouped.diagrams.len());
        for (output_index, (comparison_master, group)) in grouped
            .diagrams
            .into_iter()
            .zip(groups.iter_mut())
            .enumerate()
        {
            if group.master != output_index {
                return Err(GenerationError::InvalidGroups(format!(
                    "group {output_index} points to unexpected master {}",
                    group.master
                )));
            }
            let source_diagram = group
                .members
                .first()
                .ok_or_else(|| {
                    GenerationError::InvalidGroups(format!(
                        "group {output_index} has no representative member"
                    ))
                })?
                .source_diagram;
            let (source_representative, reversed_edges) =
                representatives.get(source_diagram).ok_or_else(|| {
                    GenerationError::InvalidGroups(format!(
                        "group {output_index} references unknown source diagram {source_diagram}"
                    ))
                })?;
            let cuts = transport_cut_metadata(comparison_master.cuts(), reversed_edges);
            let thresholds = transport_threshold_candidate_metadata(
                comparison_master.topology_threshold_candidates(),
                reversed_edges,
            );
            let mut representative = source_representative.clone();
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
                != comparison_master.loop_momentum_basis().loop_edges
            {
                representative = representative.with_loop_momentum_edges(
                    &comparison_master.loop_momentum_basis().loop_edges,
                )?;
            }
            representative = representative
                .with_name(comparison_master.name())
                .with_overall_factor(comparison_master.overall_factor().clone());
            diagrams.push(representative);

            for member in &mut group.members {
                let (source, _) = representatives.get(member.source_diagram).ok_or_else(|| {
                    GenerationError::InvalidGroups(format!(
                        "group {output_index} references unknown source diagram {}",
                        member.source_diagram
                    ))
                })?;
                member.source_id = source.id();
                member.source_name = source.name().to_owned();
            }
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

    fn validate_options(
        &self,
        generation_type: GenerationType,
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
                if options.number_of_numerical_samples == 0 && !options.check_canonical_numerator
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
                        particle: edge.data.particle,
                        direction,
                    });
                    // A self-loop contributes both half-edges to the
                    // interaction signature.
                    if edge.vertices.0 == edge.vertices.1 {
                        signature.push(EdgeColor {
                            particle: edge.data.particle,
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
                            .map(|edge| {
                                format!("{:?}:particle#{}", edge.direction, edge.particle.index())
                            })
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
                GenerationFilter::LoopOneParticleIrreducible if graph.num_loops() > 0 => {
                    let mut tree = graph.get_spanning_tree(0);
                    tree.chain_decomposition();
                    if tree.count_bridges() != 0 {
                        return Ok(false);
                    }
                }
                GenerationFilter::ParticleVeto(vetoed)
                    if graph.edges().iter().any(|edge| {
                        vetoed.iter().any(|selector| match selector {
                            ParticleSelector::Id { particle, .. } => {
                                *particle == edge.data.particle
                                    || self.model.particle_by_id(*particle).is_ok_and(|selected| {
                                        selected.antiparticle == edge.data.particle
                                    })
                            }
                            ParticleSelector::Name(_) | ParticleSelector::Pdg(_) => {
                                unreachable!("generation selectors are resolved before use")
                            }
                        })
                    }) =>
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
                | GenerationFilter::LoopOneParticleIrreducible
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
                external_particles.insert(
                    external.index,
                    self.model.particle_by_id(edge.data.particle)?,
                );
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
                let massless = self
                    .model
                    .particle_is_massless(self.model.particle_id(&particle.name)?);
                if (!massless && options.veto_massive) || (massless && options.veto_massless) {
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
                !self.model.particle_is_massless(edge.data.particle)
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
                GenerationFilter::VertexAllow(selectors) => {
                    let allowed = selectors
                        .iter()
                        .map(|selector| match selector {
                            VertexSelector::Id { vertex, .. } => *vertex,
                            VertexSelector::Name(_) => {
                                unreachable!("generation selectors are resolved before use")
                            }
                        })
                        .collect::<BTreeSet<_>>();
                    if graph.nodes().iter().any(|node| {
                        matches!(
                            node.data,
                            ColoredNode::Interaction(rule) if !allowed.contains(&rule)
                        )
                    }) {
                        return Ok(false);
                    }
                }
                GenerationFilter::VertexVeto(selectors) => {
                    let vetoed = selectors
                        .iter()
                        .map(|selector| match selector {
                            VertexSelector::Id { vertex, .. } => *vertex,
                            VertexSelector::Name(_) => {
                                unreachable!("generation selectors are resolved before use")
                            }
                        })
                        .collect::<BTreeSet<_>>();
                    if graph.nodes().iter().any(|node| {
                        matches!(
                            node.data,
                            ColoredNode::Interaction(rule) if vetoed.contains(&rule)
                        )
                    }) {
                        return Ok(false);
                    }
                }
                GenerationFilter::ParticleVeto(_)
                | GenerationFilter::MaxNumberOfBridges(_)
                | GenerationFilter::LoopOneParticleIrreducible
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
            let ColoredNode::Interaction(rule_id) = &node.data else {
                continue;
            };
            let rule = self.model.vertex_rule_by_id(*rule_id)?;
            for (name, order) in rule.coupling_orders(&self.model) {
                *result.entry(name).or_insert(0) += order;
            }
        }
        Ok(result)
    }

    // Count closed, non-branching fermion components, optionally including ghosts
    // for the Grassmann sign. The fermion-loop filter excludes ghosts. Branching
    // components are rejected explicitly below.
    fn closed_fermion_loop_count(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
        include_ghosts: bool,
    ) -> Result<usize, GenerationError> {
        let mut adjacency = vec![Vec::new(); graph.nodes().len()];
        let mut degrees = vec![0; graph.nodes().len()];
        // Build adjacency using only the selected Grassmann edges. A self-loop
        // contributes degree two.
        for edge in graph.edges() {
            let particle = self.model.particle_by_id(edge.data.particle)?;
            if !(particle.is_fermion() && !particle.is_ghost()
                || include_ghosts && particle.is_ghost())
            {
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
        generation_type: GenerationType,
        topology: ColoredTopology,
        fermion_signs: FermionSigns,
        options: &GenerationOptions,
    ) -> Result<FeynmanDiagram, GenerationError> {
        // Symbolica's automorphism factor is a denominator in the overall
        // multiplier.
        let mut factors = vec![
            symbol!("feynkit_generator_factor::AutG")
                .call(topology.symmetry)
                .pow(-1),
        ];
        if topology.multiplicity != 1 {
            factors.push(
                symbol!("feynkit_generator_factor::CouplingsMultiplicity")
                    .call(topology.multiplicity),
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
        let incoming_count = topology
            .graph
            .nodes()
            .iter()
            .filter(|node| {
                matches!(
                    node.data,
                    ColoredNode::External(ExternalNode {
                        state: ExternalState::Incoming,
                        ..
                    })
                )
            })
            .count();
        let projector =
            self.external_state_projector(&topology.graph, options.projector.as_ref())?;
        let mut builder = FeynmanDiagram::builder(Arc::clone(&self.model), name)
            .canonical_loop_momentum_basis()
            .symmetry_factor(topology.symmetry)
            .overall_factor(overall_factor)
            .numerator_prefactor(NumeratorInterface::admit(
                options.numerator_prefactor.clone(),
                None,
            )?)
            .projector(projector)
            .cuts(topology.cut_partitions.clone())
            .topology_threshold_candidates(topology.topology_threshold_candidates.clone());
        let mut endpoint_slots = BTreeMap::new();
        for (index, node) in topology.graph.nodes().iter().enumerate() {
            let vertex = match &node.data {
                ColoredNode::External(external) => {
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
                    builder.add_generation_external(
                        format!("ext{}", external.index),
                        external.index,
                        external.state,
                        connection,
                    );
                    continue;
                }
                ColoredNode::Interaction(rule_id) => {
                    let rule = self.model.vertex_rule_by_id(*rule_id)?;
                    let legs = self.interaction_half_edges(&topology.graph, index, rule)?;
                    for (slot, leg) in legs.iter().enumerate() {
                        endpoint_slots.insert((leg.edge, leg.flow), VertexSlot(slot));
                    }
                    DiagramVertex::interaction(format!("v{index}"), *rule_id)
                }
            };
            builder.add_vertex(vertex);
        }
        for (edge_index, edge) in topology.graph.edges().iter().enumerate() {
            let diagram_edge = DiagramEdge::new(edge.data.particle, edge.directed);
            let source_slot = endpoint_slots
                .get(&(edge_index, Flow::Source))
                .copied()
                .unwrap_or(VertexSlot(0));
            let target_slot = endpoint_slots
                .get(&(edge_index, Flow::Sink))
                .copied()
                .unwrap_or(VertexSlot(0));
            builder.add_edge_with_slots(
                feynkit_graph::VertexId(edge.vertices.0),
                feynkit_graph::VertexId(edge.vertices.1),
                diagram_edge,
                source_slot,
                target_slot,
            )?;
        }
        let diagram = builder.build()?;
        self.construct_numerators(diagram)
    }

    fn construct_numerators(
        &self,
        diagram: FeynmanDiagram,
    ) -> Result<FeynmanDiagram, GenerationError> {
        // Instantiate directly in the finalized frame. This avoids building
        // every factor in provisional indices and then translating it again.
        let mut vertices = Vec::new();
        for (vertex, data) in diagram.vertices() {
            let Some(rule_id) = data.interaction else {
                continue;
            };
            let rule = self.model.vertex_rule_by_id(rule_id)?;
            let mut legs = vec![None; rule.particles.len()];
            for (edge, endpoints, data) in diagram.edges() {
                for (attached, slot, endpoint, flow) in [
                    (
                        endpoints.source,
                        data.source_slot(),
                        DiagramEndpoint::Source,
                        Flow::Source,
                    ),
                    (
                        endpoints.target,
                        data.target_slot(),
                        DiagramEndpoint::Target,
                        Flow::Sink,
                    ),
                ] {
                    if attached != Some(vertex) {
                        continue;
                    }
                    let particle = self.model.particle_by_id(rule.particles[slot.0])?;
                    let hedge = diagram
                        .half_edge_id(DiagramHalfEdge { edge, endpoint })
                        .expect("finalized endpoint")
                        .0;
                    legs[slot.0] = Some(NumeratorHalfEdge {
                        edge: edge.0,
                        hedge,
                        flow,
                        spin: particle.spin,
                        color: particle.color,
                    });
                }
            }
            let legs = legs
                .into_iter()
                .map(|leg| leg.expect("assigned interaction slot"))
                .collect::<Vec<_>>();
            let instantiation = NumeratorInstantiation {
                owner: NumeratorOwner::Vertex(vertex.0),
                legs: &legs,
            };
            vertices.push((vertex, self.interaction_numerator(rule_id, &instantiation)?));
        }
        let mut edges = Vec::new();
        for (edge, _, data) in diagram
            .edges()
            .filter(|(_, _, data)| data.external.is_none())
        {
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

    fn external_state_projector(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
        supplied: Option<&Atom>,
    ) -> Result<SymbolicTensor<PartialStructure>, GenerationError> {
        let mut projector = Atom::one();
        let mut external_indices = BTreeMap::new();
        for (vertex, node) in graph.nodes().iter().enumerate() {
            let ColoredNode::External(external) = &node.data else {
                continue;
            };
            let edge_id = *node
                .edges
                .first()
                .ok_or(GenerationError::MissingExternalEdge { node: vertex })?;
            let edge = &graph.edges()[edge_id];
            let flow = if edge.vertices.0 == vertex {
                Flow::Sink
            } else {
                Flow::Source
            };
            let hedge = edge_id
                .checked_mul(2)
                .and_then(|index| index.checked_add(usize::from(flow == Flow::Sink)))
                .ok_or(GenerationError::ArithmeticOverflow(
                    "an external half-edge identifier",
                ))?;
            if supplied.is_some() {
                external_indices.insert(external.index, (edge_id, hedge));
            } else {
                projector *=
                    self.external_wavefunction(external.state, external.particle, edge_id, hedge)?;
            }
        }
        let Some(supplied) = supplied else {
            return Ok(NumeratorInterface::admit(projector, None)?);
        };
        // Process-level expressions cannot know a topology's provisional IDs.
        // Resolve their external-leg owners before the finalizer transports the
        // resulting edge/half-edge chart. Native local UFO slots are unchanged.
        let hedge_head = symbols::hedge_index();
        let edge_head = symbols::edge_index();
        let carriers = [
            momentum_symbol(),
            symbols::u(),
            symbols::ubar(),
            symbols::v(),
            symbols::vbar(),
            symbols::epsilon(),
            symbols::epsilonbar(),
        ];
        projector = supplied.clone();
        for carrier_pass in [false, true] {
            let mut unknown = None;
            projector = projector.replace_map(|part, _, output| {
                let AtomView::Fun(function) = part else {
                    return;
                };
                let head = function.get_symbol();
                let is_index = head == hedge_head || head == edge_head;
                if (carrier_pass && !carriers.contains(&head)) || (!carrier_pass && !is_index) {
                    return;
                }
                let Some(index) = function
                    .iter()
                    .next()
                    .and_then(|arg| usize::try_from(arg).ok())
                else {
                    return;
                };
                let Some(&(edge, hedge)) = external_indices.get(&index) else {
                    unknown = Some(index);
                    return;
                };
                **output = FunctionBuilder::new(head)
                    .add_arg(if head == hedge_head { hedge } else { edge })
                    .add_args(function.iter().skip(1))
                    .finish();
            });
            if let Some(index) = unknown {
                return Err(GenerationError::UnknownProjectorExternalIndex { index });
            }
        }
        Ok(NumeratorInterface::admit(projector, None)?)
    }

    fn external_wavefunction(
        &self,
        state: ExternalState,
        particle: ParticleId,
        edge: usize,
        hedge: usize,
    ) -> Result<Atom, GenerationError> {
        let index = symbols::hedge_index().call((hedge, 1_i64));
        let particle = self.model.particle_by_id(particle)?;
        let (head, index) = match particle.spin {
            2 => {
                let head = match (state, particle.is_antiparticle()) {
                    (ExternalState::Incoming, false) => symbols::u(),
                    (ExternalState::Incoming, true) => symbols::vbar(),
                    (ExternalState::Outgoing, false) => symbols::ubar(),
                    (ExternalState::Outgoing, true) => symbols::v(),
                };
                (head, Bispinor {}.new_rep(4).to_symbolic([index]))
            }
            3 => {
                let head = match state {
                    ExternalState::Incoming => symbols::epsilon(),
                    ExternalState::Outgoing => symbols::epsilonbar(),
                };
                (head, Minkowski {}.new_rep(4).to_symbolic([index]))
            }
            _ => return Ok(Atom::one()),
        };
        Ok(FunctionBuilder::new(head)
            .add_arg(edge)
            .add_arg(index)
            .finish())
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
        let mut available = available_interaction_half_edges(graph, vertex, &self.model)?;

        let mut ordered = Vec::with_capacity(rule.particles.len());
        for (leg, particle_id) in rule.particles.iter().enumerate() {
            let particle = self.model.particle_by_id(*particle_id)?;
            let expected = ResolvedProcess::vertex_half_edge(&self.model, particle)?;
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
    ) -> Result<SymbolicTensor<PartialStructure>, GenerationError> {
        let key = InteractionNumeratorKey {
            rule: rule_id,
            owner: instantiation.owner,
            legs: instantiation.legs.to_vec(),
        };
        if let Some(numerator) = self
            .numerators
            .interactions
            .read()
            .expect("numerator cache was poisoned; its retained tensors cannot be reused")
            .get(&key)
        {
            return Ok(numerator.clone());
        }
        let legs = instantiation
            .legs
            .iter()
            .enumerate()
            .map(|(index, leg)| NumeratorHalfEdge {
                edge: instantiation.legs[..index]
                    .iter()
                    .position(|previous| previous.edge == leg.edge)
                    .unwrap_or(index),
                hedge: index,
                ..*leg
            })
            .collect::<Vec<_>>();
        let template = self.interaction_template(rule_id, legs)?;
        let numerator = match &*template {
            NumeratorTemplate::Dynamic => NumeratorInterface::admit(
                self.interaction_expression(rule_id, instantiation)?,
                None,
            )?,
            NumeratorTemplate::Admitted { .. } => {
                let numerator = template.instantiate(instantiation)?;
                // Raw model rules have no declared axis order: retain the
                // normalized occurrence order used by actual-owner admission.
                let structure = composition::interface_in_atom_order(
                    numerator.expression(),
                    numerator.structure(),
                )
                .expect("template admission established enumerable explicit ports");
                let logical = numerator.structure().logical_slots();
                let axes = structure
                    .logical_slots()
                    .iter()
                    .map(|slot| {
                        logical
                            .iter()
                            .position(|original| original == slot)
                            .expect("template transport preserves every explicit interface slot")
                    })
                    .collect::<Vec<_>>();
                numerator.permuted(&axes).map_err(DiagramError::from)?
            }
        };
        self.numerators
            .interactions
            .write()
            .expect("numerator cache was poisoned; its retained tensors cannot be reused")
            .insert(key, numerator.clone());
        Ok(numerator)
    }

    fn interaction_template(
        &self,
        rule_id: VertexRuleId,
        legs: Vec<NumeratorHalfEdge>,
    ) -> Result<Arc<NumeratorTemplate>, GenerationError> {
        let key = InteractionNumeratorKey {
            rule: rule_id,
            owner: NumeratorOwner::Vertex(0),
            legs,
        };
        if let Some(template) = self
            .numerators
            .interaction_templates
            .read()
            .expect("numerator cache was poisoned; its retained tensors cannot be reused")
            .get(&key)
        {
            return Ok(template.clone());
        }
        let distinct = key
            .legs
            .iter()
            .enumerate()
            .map(|(index, leg)| NumeratorHalfEdge {
                edge: index,
                ..*leg
            })
            .collect::<Vec<_>>();
        // Establish the unspecialized rule first. Aliasing momenta can cancel a
        // tensor to zero, but must not hide malformed summands or erase its rank.
        let base = if distinct != key.legs {
            Some(self.interaction_template(rule_id, distinct)?)
        } else {
            None
        };
        let instantiation = NumeratorInstantiation {
            owner: key.owner,
            legs: &key.legs,
        };
        let template = if let Some(base) = base {
            match &*base {
                NumeratorTemplate::Dynamic => NumeratorTemplate::Dynamic,
                NumeratorTemplate::Admitted { tensor, labels } => {
                    let expression = instantiation.bind_template(tensor.expression())?;
                    let specialized = tensor
                        .with_renamed_indices(expression, &HashMap::new())
                        .map_err(DiagramError::from)?;
                    // Momentum aliases can erase the payload. Keep the admitted
                    // rule's complete label universe for its retained metadata.
                    NumeratorTemplate::Admitted {
                        tensor: Box::new(specialized),
                        labels: labels.clone(),
                    }
                }
            }
        } else {
            let rule = self.model.vertex_rule_by_id(rule_id)?;
            let mut supported = rule
                .color_structures
                .iter()
                .all(NumeratorTemplate::supports);
            for lorentz in &rule.lorentz_structures {
                supported &= NumeratorTemplate::supports(
                    &self.model.lorentz_structure_by_id(*lorentz)?.structure,
                );
            }
            if supported {
                let expression = self.interaction_expression(rule_id, &instantiation)?;
                NumeratorTemplate::admitted(NumeratorInterface::admit(expression, None)?)?
            } else {
                NumeratorTemplate::Dynamic
            }
        };
        // Build outside the lock: symbolic normalization may invoke callbacks.
        // Concurrent cold callers may compute the same immutable rule; only the
        // first completed template is retained and returned to all callers.
        let mut templates = self
            .numerators
            .interaction_templates
            .write()
            .expect("numerator template cache was poisoned while publishing a rule");
        Ok(Arc::clone(
            templates.entry(key).or_insert_with(|| Arc::new(template)),
        ))
    }

    fn interaction_expression(
        &self,
        rule_id: VertexRuleId,
        instantiation: &NumeratorInstantiation<'_>,
    ) -> Result<Atom, GenerationError> {
        let rule = self.model.vertex_rule_by_id(rule_id)?;
        let mut terms = Vec::new();
        for (color_index, color) in rule.color_structures.iter().enumerate() {
            for (lorentz_index, lorentz_name) in rule.lorentz_structures.iter().enumerate() {
                if let Some(coupling) = rule
                    .couplings
                    .get(color_index)
                    .and_then(|row| row.get(lorentz_index))
                    .and_then(Option::as_ref)
                {
                    let coupling = self.model.coupling_by_id(*coupling)?;
                    let coupling = Atom::var(symbol!(&format!("UFO::{}", coupling.name)));
                    let lorentz = &self.model.lorentz_structure_by_id(*lorentz_name)?.structure;
                    let color = instantiation.instantiate(color, NumeratorSector::Color)?;
                    let lorentz = instantiation.instantiate(lorentz, NumeratorSector::Spin)?;
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

    fn edge_numerator(
        &self,
        particle_id: ParticleId,
        edge: usize,
        hedges: [usize; 2],
    ) -> Result<SymbolicTensor<PartialStructure>, GenerationError> {
        let key = PropagatorNumeratorKey {
            particle: particle_id,
            edge,
            hedges,
        };
        if let Some(numerator) = self
            .numerators
            .propagators
            .read()
            .expect("numerator cache was poisoned; its retained tensors cannot be reused")
            .get(&key)
        {
            return Ok(numerator.clone());
        }
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
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Edge(edge),
            legs: &legs,
        };
        let cached = self
            .numerators
            .propagator_templates
            .read()
            .expect("numerator cache was poisoned; its retained tensors cannot be reused")
            .get(&particle_id)
            .cloned();
        let template = if let Some(template) = cached {
            template
        } else {
            let template = if NumeratorTemplate::supports(&propagator) {
                let canonical = legs
                    .iter()
                    .enumerate()
                    .map(|(index, leg)| NumeratorHalfEdge {
                        edge: 0,
                        hedge: index,
                        ..*leg
                    })
                    .collect::<Vec<_>>();
                let instantiation = NumeratorInstantiation {
                    owner: NumeratorOwner::Edge(0),
                    legs: &canonical,
                };
                let expression =
                    instantiation.propagator_expression(&propagator, particle.color)?;
                NumeratorTemplate::admitted(NumeratorInterface::admit(expression, None)?)?
            } else {
                NumeratorTemplate::Dynamic
            };
            // Publish only after normalization, without holding a lock across
            // a callback. Reuse the winner if another cold caller raced us.
            let mut templates = self
                .numerators
                .propagator_templates
                .write()
                .expect("propagator template cache was poisoned while publishing a rule");
            Arc::clone(
                templates
                    .entry(particle_id)
                    .or_insert_with(|| Arc::new(template)),
            )
        };
        let numerator = match &*template {
            NumeratorTemplate::Dynamic => NumeratorInterface::admit(
                instantiation.propagator_expression(&propagator, particle.color)?,
                None,
            )?,
            NumeratorTemplate::Admitted { .. } => {
                let numerator = template.instantiate(&instantiation)?;
                let structure = composition::interface_in_atom_order(
                    numerator.expression(),
                    numerator.structure(),
                )
                .expect("template admission established enumerable explicit ports");
                let logical = numerator.structure().logical_slots();
                let axes = structure
                    .logical_slots()
                    .iter()
                    .map(|slot| {
                        logical
                            .iter()
                            .position(|original| original == slot)
                            .expect("template transport preserves every explicit interface slot")
                    })
                    .collect::<Vec<_>>();
                numerator.permuted(&axes).map_err(DiagramError::from)?
            }
        };
        self.numerators
            .propagators
            .write()
            .expect("numerator cache was poisoned; its retained tensors cannot be reused")
            .insert(key, numerator.clone());
        Ok(numerator)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdgeColor {
    pub particle: ParticleId,
    pub direction: Option<bool>,
}

impl fmt::Display for EdgeColor {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "particle#{}", self.particle.index())
    }
}

fn available_interaction_half_edges(
    graph: &Graph<ColoredNode, EdgeColor>,
    vertex: usize,
    model: &Model,
) -> Result<Vec<(NumeratorHalfEdge, HalfEdge<EdgeColor>)>, GenerationError> {
    let mut available = Vec::new();
    for edge_id in &graph.nodes()[vertex].edges {
        let edge = &graph.edges()[*edge_id];
        let particle = model.particle_by_id(edge.data.particle)?;
        let base = if particle.is_antiparticle() {
            particle.antiparticle
        } else {
            edge.data.particle
        };
        if edge.vertices.0 == vertex {
            available.push((
                NumeratorHalfEdge {
                    edge: *edge_id,
                    flow: Flow::Source,
                    hedge: 2 * *edge_id,
                    spin: 0,
                    color: 1,
                },
                HalfEdge {
                    data: EdgeColor {
                        particle: base,
                        direction: None,
                    },
                    direction: edge.directed.then_some(!particle.is_antiparticle()),
                },
            ));
        }
        if edge.vertices.1 == vertex {
            available.push((
                NumeratorHalfEdge {
                    edge: *edge_id,
                    flow: Flow::Sink,
                    hedge: 2 * *edge_id + 1,
                    spin: 0,
                    color: 1,
                },
                HalfEdge {
                    data: EdgeColor {
                        particle: base,
                        direction: None,
                    },
                    direction: edge.directed.then_some(particle.is_antiparticle()),
                },
            ));
        }
    }
    Ok(available)
}

fn vertex_rule_signature(
    model: &Model,
    rule: &VertexRule,
) -> Result<Vec<HalfEdge<EdgeColor>>, GenerationError> {
    let mut signature = rule
        .particles
        .iter()
        .map(|particle| ResolvedProcess::vertex_half_edge(model, model.particle_by_id(*particle)?))
        .collect::<Result<Vec<_>, GenerationError>>()?;
    signature.sort();
    Ok(signature)
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ExternalNode {
    pub index: usize,
    pub state: ExternalState,
    pub particle: ParticleId,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum ExternalCanonicalClass {
    Exact { state: ExternalState, index: usize },
    WithinState(ExternalState),
}

#[derive(Debug, Clone, Default)]
pub struct NodeColor {
    pub external: Option<ExternalNode>,
    external_class: Option<ExternalCanonicalClass>,
}

impl NodeColor {
    fn key(&self) -> Option<(ParticleId, ExternalCanonicalClass)> {
        self.external.as_ref().map(|external| {
            (
                external.particle,
                self.external_class
                    .expect("external nodes always have a canonical class"),
            )
        })
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
    Interaction(VertexRuleId),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct NumeratorCanonicalNode(i64);

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct NumeratorCanonicalEdge(String);

type InteractionAssignment = (Graph<ColoredNode, EdgeColor>, usize);

impl fmt::Display for ColoredNode {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::External(external) => write!(formatter, "ext{}", external.index),
            Self::Interaction(rule) => write!(formatter, "vertex#{}", rule.index()),
        }
    }
}

#[derive(Clone, PartialEq, Eq)]
struct ColoredTopology {
    graph: Graph<ColoredNode, EdgeColor>,
    symmetry: u64,
    multiplicity: usize,
    cut_partitions: Vec<DiagramCut>,
    topology_threshold_candidates: Vec<DiagramThresholdCandidate>,
}

struct DiagramPair {
    comparison: FeynmanDiagram,
    representative: FeynmanDiagram,
    reversed_edges: BTreeSet<EdgeId>,
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

#[derive(Debug, Clone, Copy, Default)]
struct FermionSigns {
    internal_loop_negative: bool,
    include_external_ordering: bool,
    external_ordering_negative: bool,
    include_antifermion_spin_sum: bool,
    antifermion_spin_sum_negative: bool,
}

struct NormalizedFermionFlow {
    #[cfg_attr(not(test), allow(dead_code))]
    graph: Graph<ColoredNode, EdgeColor>,
    reversed_edges: BTreeSet<EdgeId>,
    signs: FermionSigns,
}

struct FermionNormalizationState<'a> {
    visited_edges: &'a mut [bool],
    vertices: &'a mut [Option<(usize, usize)>],
    particles: &'a mut [Option<ParticleId>],
}

fn normalize_fermion_component(
    graph: &Graph<ColoredNode, EdgeColor>,
    starting_edge: usize,
    adjacency: &[Vec<(usize, usize)>],
    state: &mut FermionNormalizationState<'_>,
    model: &Model,
    canonicalize_virtual_start: bool,
) -> Result<BTreeSet<usize>, GenerationError> {
    let starting = &graph.edges()[starting_edge];
    let starting_particle = model.particle_by_id(starting.data.particle)?;
    let reverse_start = canonicalize_virtual_start && starting_particle.is_antiparticle();
    let (left, right) = if reverse_start {
        (starting.vertices.1, starting.vertices.0)
    } else {
        starting.vertices
    };
    let chain_is_antiparticle = starting_particle.is_antiparticle() && !reverse_start;
    state.visited_edges[starting_edge] = true;
    state.vertices[starting_edge] = Some((left, right));
    state.particles[starting_edge] = Some(if reverse_start {
        starting_particle.antiparticle
    } else {
        starting.data.particle
    });
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
            .find(|(edge_id, _)| !state.visited_edges[*edge_id])
        {
            state.visited_edges[edge_id] = true;
            let vertices = if read_forward {
                (current, next)
            } else {
                (next, current)
            };
            state.vertices[edge_id] = Some(vertices);
            let edge = &graph.edges()[edge_id];
            let particle = model.particle_by_id(edge.data.particle)?;
            // Legacy GammaLoop represents a complete charged chain with the
            // same particle (or antiparticle) on every edge.  When a
            // canonical edge points against the chain, reversing only its
            // endpoints is not enough: choose the conjugate precisely when
            // needed to match the species selected by the starting edge.
            state.particles[edge_id] = Some(
                if vertices != edge.vertices && particle.is_antiparticle() != chain_is_antiparticle
                {
                    particle.antiparticle
                } else {
                    edge.data.particle
                },
            );
            record_external(next);
            current = next;
        }
    }
    Ok(external_legs)
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

struct InteractionSignatures(BTreeMap<Vec<HalfEdge<EdgeColor>>, Vec<VertexRuleId>>);

impl InteractionSignatures {
    fn new(model: &Model, options: &GenerationOptions) -> Result<Self, GenerationError> {
        let allowed = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::VertexAllow(selectors) => Some(
                    selectors
                        .iter()
                        .map(|selector| match selector {
                            VertexSelector::Id { vertex, .. } => *vertex,
                            VertexSelector::Name(_) => {
                                unreachable!("generation selectors are resolved before use")
                            }
                        })
                        .collect::<BTreeSet<_>>(),
                ),
                _ => None,
            });
        let vetoed = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::VertexVeto(selectors) => Some(
                    selectors
                        .iter()
                        .map(|selector| match selector {
                            VertexSelector::Id { vertex, .. } => *vertex,
                            VertexSelector::Name(_) => {
                                unreachable!("generation selectors are resolved before use")
                            }
                        })
                        .collect::<BTreeSet<_>>(),
                ),
                _ => None,
            });
        let particle_vetoes = options
            .graph_filters
            .iter()
            .find_map(|filter| match filter {
                GenerationFilter::ParticleVeto(selectors) => Some(
                    selectors
                        .iter()
                        .map(|selector| match selector {
                            ParticleSelector::Id { particle, .. } => *particle,
                            ParticleSelector::Name(_) | ParticleSelector::Pdg(_) => {
                                unreachable!("generation selectors are resolved before use")
                            }
                        })
                        .collect::<BTreeSet<_>>(),
                ),
                _ => None,
            });
        let mut signatures: BTreeMap<Vec<_>, Vec<_>> = BTreeMap::new();
        for (rule_index, rule) in model.vertex_rules().iter().enumerate() {
            let rule_id = model.vertex_rule_id_at(rule_index)?;
            if allowed
                .as_ref()
                .is_some_and(|allowed| !allowed.contains(&rule_id))
                || vetoed
                    .as_ref()
                    .is_some_and(|vetoed| vetoed.contains(&rule_id))
            {
                continue;
            }
            // Coupling powers are nonnegative. A vertex exceeding a graph-wide
            // upper bound cannot occur in an accepted graph, so exclude its
            // signature before enumerating topologies (especially QED=0).
            let orders = rule.coupling_orders(model);
            if options.graph_filters.iter().any(|filter| {
                matches!(filter, GenerationFilter::CouplingOrders(required)
                    if required.iter().any(|(name, (_, maximum))|
                        maximum.is_some_and(|maximum| orders.get(name).copied().unwrap_or(0) > maximum)))
            }) {
                continue;
            }
            let mut signature = Vec::new();
            let mut rejected = false;
            for particle_id in &rule.particles {
                let particle = model.particle_by_id(*particle_id)?;
                if particle_vetoes.as_ref().is_some_and(|vetoed| {
                    vetoed.contains(particle_id) || vetoed.contains(&particle.antiparticle)
                }) {
                    rejected = true;
                    break;
                }
                signature.push(ResolvedProcess::vertex_half_edge(model, particle)?);
            }
            if !rejected {
                signature.sort();
                signatures.entry(signature).or_default().push(rule_id);
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

    fn rules(&self, signature: &[EdgeColor]) -> Option<&[VertexRuleId]> {
        let mut half_edges: Vec<_> = signature
            .iter()
            .map(|edge| HalfEdge {
                direction: edge.direction,
                data: EdgeColor {
                    particle: edge.particle,
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
            particle.antiparticle
        } else {
            model.particle_id(&particle.name)?
        };
        let particle_id = model.particle_id(&particle.name)?;
        let direction = if model.particle_is_self_conjugate(particle_id) {
            None
        } else {
            Some(match (state, particle.is_antiparticle()) {
                (ExternalState::Incoming, false) | (ExternalState::Outgoing, true) => true,
                (ExternalState::Incoming, true) | (ExternalState::Outgoing, false) => false,
            })
        };
        let symmetrized_within_state = self.generation_type == GenerationType::Amplitude
            && match state {
                ExternalState::Incoming => self.symmetrize_initial,
                ExternalState::Outgoing => self.symmetrize_final,
            }
            && (!particle.is_fermion() || self.symmetrize_external_fermions);
        // Both sides of a forward graph are labelled copies of the incoming
        // state. Neither independent permutations nor final-state symmetry may
        // enter its raw automorphism denominator. Paired incoming permutations
        // and the optional global side exchange are compared after generation.
        let external_class = if symmetrized_within_state {
            ExternalCanonicalClass::WithinState(state)
        } else {
            ExternalCanonicalClass::Exact { state, index }
        };
        Ok((
            NodeColor {
                external: Some(ExternalNode {
                    index,
                    state,
                    particle: model.particle_id(&particle.name)?,
                }),
                external_class: Some(external_class),
            },
            HalfEdge {
                direction,
                data: EdgeColor {
                    particle: base,
                    direction: None,
                },
            },
        ))
    }

    /// Choose the reduced-skeleton representative before any numerator
    /// indices are instantiated.
    ///
    /// For left/right symmetry the extra orbit operation is global: either
    /// every sewn incoming leg is exchanged with its outgoing partner, or none
    /// is. Keeping the exact external colors in each skeleton prevents graph
    /// automorphisms from swapping only a subset of the pairs. Once the
    /// smallest skeleton is selected, the complete colored graph is rebuilt in
    /// that canonical vertex/edge order so UFO slot matching and propagator
    /// momenta are instantiated identically for isomorphic diagrams.
    fn canonicalize_numerator_graph(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
        model: &Model,
        grouping: &NumeratorGrouping,
    ) -> Result<Graph<ColoredNode, EdgeColor>, GenerationError> {
        let grouping_options = match grouping {
            NumeratorGrouping::Identical(options)
            | NumeratorGrouping::UpToSign(options)
            | NumeratorGrouping::UpToScalar(options) => Some(options),
            NumeratorGrouping::None => None,
        };
        let include_external = self.generation_type == GenerationType::CrossSection
            && (self.symmetrize_initial || self.symmetrize_left_right);
        let exact_internal_particles =
            grouping_options.is_some_and(|options| !options.differentiate_particle_masses_only);

        let selected = self.canonicalize_numerator_graph_pass(
            graph,
            model,
            exact_internal_particles,
            include_external,
            include_external,
        )?;
        let canonical = if include_external {
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
        graph: &Graph<ColoredNode, EdgeColor>,
        model: &Model,
        exact_internal_particles: bool,
        include_external_remapping: bool,
        use_sewn_external_colors: bool,
    ) -> Result<Graph<ColoredNode, EdgeColor>, GenerationError> {
        let mut candidates = Vec::new();
        let mut remappings = vec![(false, BTreeMap::new())];
        if include_external_remapping {
            let external_nodes = graph
                .nodes()
                .iter()
                .enumerate()
                .filter_map(|(node, value)| match &value.data {
                    ColoredNode::External(external) => Some((external.index, node)),
                    ColoredNode::Interaction(_) => None,
                })
                .collect::<BTreeMap<_, _>>();
            let incoming_count = self.incoming.len();
            if external_nodes.len() != 2 * incoming_count {
                return Err(GenerationError::GraphConstruction(format!(
                    "sewn canonicalization expected {} external nodes, found {}",
                    2 * incoming_count,
                    external_nodes.len()
                )));
            }
            let mut generators = vec![Permutation::id(incoming_count)];
            if self.symmetrize_initial {
                let mut species = BTreeMap::<i64, Vec<usize>>::new();
                for (index, particle) in self.incoming.iter().enumerate() {
                    species.entry(particle.pdg_code).or_default().push(index);
                }
                // Adjacent transpositions generate only permutations within
                // each identical-species class, without enumerating rejected
                // permutations of distinct incoming particles.
                for positions in species.values() {
                    for pair in positions.windows(2) {
                        let mut mapping = (0..incoming_count).collect::<Vec<_>>();
                        mapping.swap(pair[0], pair[1]);
                        generators.push(Permutation::from_map(mapping));
                    }
                }
            }
            let mut permutations = Permutation::generate_all(&generators)
                .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
            permutations.sort_by(|left, right| left.map().cmp(right.map()));
            remappings.clear();
            for permutation in permutations {
                for swapped in [false, true] {
                    if swapped && !self.symmetrize_left_right {
                        continue;
                    }
                    let mapping = (0..incoming_count)
                        .flat_map(|index| {
                            let target = permutation.map()[index];
                            let (left, right) = if swapped {
                                (target + incoming_count, target)
                            } else {
                                (target, target + incoming_count)
                            };
                            [
                                (external_nodes[&index], external_nodes[&left]),
                                (
                                    external_nodes[&(index + incoming_count)],
                                    external_nodes[&right],
                                ),
                            ]
                        })
                        .collect();
                    remappings.push((swapped, mapping));
                }
            }
        }
        for (swapped, remapping) in remappings {
            let mut skeleton = Graph::new();
            for node in graph.nodes() {
                let color = match &node.data {
                    ColoredNode::External(external) => {
                        let index = i64::try_from(external.index).map_err(|_| {
                            GenerationError::ArithmeticOverflow("an external canonical color")
                        })?;
                        let incoming_count = i64::try_from(self.incoming.len()).map_err(|_| {
                            GenerationError::ArithmeticOverflow("an external canonical color")
                        })?;
                        if use_sewn_external_colors {
                            if index < incoming_count {
                                -2001 - index
                            } else {
                                -3001 - (index - incoming_count)
                            }
                        } else {
                            index + 1
                        }
                    }
                    ColoredNode::Interaction(_) => 0,
                };
                skeleton.add_node(NumeratorCanonicalNode(color));
            }
            for edge in graph.edges() {
                let external = matches!(
                    graph.nodes()[edge.vertices.0].data,
                    ColoredNode::External(_)
                ) || matches!(
                    graph.nodes()[edge.vertices.1].data,
                    ColoredNode::External(_)
                );
                let mut vertices = (
                    remapping
                        .get(&edge.vertices.0)
                        .copied()
                        .unwrap_or(edge.vertices.0),
                    remapping
                        .get(&edge.vertices.1)
                        .copied()
                        .unwrap_or(edge.vertices.1),
                );
                // As in the legacy comparison graph, all external edges point
                // into the interaction graph. Internal fermion direction is
                // omitted from the reduced skeleton but restored below.
                let external_outgoing = matches!(
                    graph.nodes()[edge.vertices.1].data,
                    ColoredNode::External(_)
                );
                if external_outgoing {
                    vertices = (vertices.1, vertices.0);
                }
                let particle = model.particle_by_id(if external_outgoing {
                    model.particle_by_id(edge.data.particle)?.antiparticle
                } else {
                    edge.data.particle
                })?;
                // Ghost vertices differentiate one propagator end. Preserve its
                // physical arrow in the momentum chart: reversing storage while
                // conjugating the ghost must describe the same directed edge.
                if particle.is_ghost() && particle.is_antiparticle() {
                    vertices = (vertices.1, vertices.0);
                }
                let keep_direction = external || particle.is_ghost();
                let color = if external || exact_internal_particles {
                    particle.name.clone()
                } else {
                    format!(
                        "{} | {}",
                        model.parameter_by_id(particle.mass)?.name,
                        particle.spin
                    )
                };
                skeleton
                    .add_edge(
                        vertices.0,
                        vertices.1,
                        keep_direction,
                        NumeratorCanonicalEdge(color),
                    )
                    .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
            }
            candidates.push((swapped, remapping, skeleton.canonize()));
        }
        candidates.sort_by(|left, right| {
            (left.2.graph.nodes(), left.2.graph.edges())
                .cmp(&(right.2.graph.nodes(), right.2.graph.edges()))
        });
        let (left_right_swapped, remapping, canonical) = candidates
            .into_iter()
            .next()
            .ok_or_else(|| GenerationError::GraphConstruction("no canonical graph".to_owned()))?;
        let inverse_remapping = remapping
            .iter()
            .map(|(source, target)| (*target, *source))
            .collect::<BTreeMap<_, _>>();
        let mut canonical_to_input = vec![0; graph.nodes().len()];
        for (input, canonical_position) in canonical.vertex_map.iter().copied().enumerate() {
            canonical_to_input[canonical_position] =
                inverse_remapping.get(&input).copied().unwrap_or(input);
        }
        let mut input_to_canonical = vec![0; graph.nodes().len()];
        for (canonical_position, input) in canonical_to_input.iter().copied().enumerate() {
            input_to_canonical[input] = canonical_position;
        }

        let mut result = Graph::new();
        for input in canonical_to_input.iter().copied() {
            // Source edges retain their input owners; external labels belong
            // to the selected paired remapping. General incoming cycles need
            // the inverse above, unlike the involutive global side exchange.
            let data = remapping.get(&input).copied().unwrap_or(input);
            result.add_node(graph.nodes()[data].data.clone());
        }

        let mut edges = graph
            .edges()
            .iter()
            .map(|edge| {
                let external = matches!(
                    graph.nodes()[edge.vertices.0].data,
                    ColoredNode::External(_)
                ) || matches!(
                    graph.nodes()[edge.vertices.1].data,
                    ColoredNode::External(_)
                );
                let external_outgoing = matches!(
                    graph.nodes()[edge.vertices.1].data,
                    ColoredNode::External(_)
                );
                let flipped = external_outgoing
                    || (!external
                        && input_to_canonical[edge.vertices.0]
                            > input_to_canonical[edge.vertices.1]);
                let vertices = if flipped {
                    (
                        input_to_canonical[edge.vertices.1],
                        input_to_canonical[edge.vertices.0],
                    )
                } else {
                    (
                        input_to_canonical[edge.vertices.0],
                        input_to_canonical[edge.vertices.1],
                    )
                };
                let mut data = edge.data;
                if edge.directed && flipped {
                    data.particle = model.particle_by_id(data.particle)?.antiparticle;
                }
                let particle = model.particle_by_id(data.particle)?;
                // CP acts on neutral adjoint Faddeev-Popov ghosts through
                // color, without exchanging ghost and antighost fields. Their
                // derivative momentum must keep its ghost-number orientation.
                if edge.directed
                    && left_right_swapped
                    && !(particle.is_ghost() && particle.color == 8 && particle.charge.is_zero())
                {
                    data.particle = particle.antiparticle;
                }
                let pdg_code = model.particle_by_id(data.particle)?.pdg_code;
                Ok((vertices, edge.directed, data, pdg_code))
            })
            .collect::<Result<Vec<_>, GenerationError>>()?;
        // Match the UFO/GammaLoop canonical momentum convention. Signed PDG
        // order is deliberately used here rather than model insertion order:
        // CP-conjugate charged lines can otherwise acquire the same edge/LMB
        // order even when their legacy canonical representatives differ.
        edges.sort_by_key(|(vertices, _, _, pdg_code)| (vertices.0, vertices.1, *pdg_code));
        for (vertices, directed, data, _) in edges {
            result
                .add_edge(vertices.0, vertices.1, directed, data)
                .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
        }
        Ok(result)
    }

    fn repair_cp_interaction_rules(
        graph: Graph<ColoredNode, EdgeColor>,
        model: &Model,
    ) -> Result<Graph<ColoredNode, EdgeColor>, GenerationError> {
        let mut replacements = BTreeMap::new();
        for (vertex, node) in graph.nodes().iter().enumerate() {
            let ColoredNode::Interaction(rule_id) = &node.data else {
                continue;
            };
            let rule = model.vertex_rule_by_id(*rule_id)?;
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
                .map(|half_edge| HalfEdge {
                    direction: half_edge.direction.map(|direction| !direction),
                    data: half_edge.data,
                })
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
        let mut repaired = Graph::new();
        for (vertex, node) in graph.nodes().iter().enumerate() {
            let data = replacements
                .get(&vertex)
                .map_or_else(|| node.data.clone(), |rule| ColoredNode::Interaction(*rule));
            repaired.add_node(data);
        }
        for edge in graph.edges() {
            repaired
                .add_edge(edge.vertices.0, edge.vertices.1, edge.directed, edge.data)
                .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
        }
        Ok(repaired)
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
            let particle = model.particle_by_id(edge.data.particle)?;
            if !(particle.is_fermion() || particle.is_ghost()) {
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
        let mut normalized_particles = vec![None; graph.edges().len()];
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
            let Some((edge_id, _)) = adjacency[*node_id]
                .iter()
                .find(|(edge_id, _)| !normalization.visited_edges[*edge_id])
            else {
                continue;
            };
            let legs = normalize_fermion_component(
                graph,
                *edge_id,
                &adjacency,
                &mut normalization,
                model,
                false,
            )?;
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
        for edge_id in 0..graph.edges().len() {
            let particle = model.particle_by_id(graph.edges()[edge_id].data.particle)?;
            if (particle.is_fermion() || particle.is_ghost())
                && !normalization.visited_edges[edge_id]
            {
                normalize_fermion_component(
                    graph,
                    edge_id,
                    &adjacency,
                    &mut normalization,
                    model,
                    true,
                )?;
            }
        }

        let mut normalized_graph = Graph::new();
        for node in graph.nodes() {
            normalized_graph.add_node(node.data.clone());
        }
        let mut reversed_edges = BTreeSet::new();
        for (edge_id, edge) in graph.edges().iter().enumerate() {
            let particle = model.particle_by_id(edge.data.particle)?;
            let canonicalize_virtual_antiparticle = normalized_vertices[edge_id].is_none()
                && matches!(
                    (
                        &graph.nodes()[edge.vertices.0].data,
                        &graph.nodes()[edge.vertices.1].data
                    ),
                    (ColoredNode::Interaction(_), ColoredNode::Interaction(_))
                )
                && particle.is_antiparticle()
                && !(particle.is_fermion() || particle.is_ghost());
            let (left, right) = normalized_vertices[edge_id].unwrap_or({
                if canonicalize_virtual_antiparticle {
                    (edge.vertices.1, edge.vertices.0)
                } else {
                    edge.vertices
                }
            });
            let mut data = edge.data;
            data.particle = normalized_particles[edge_id].unwrap_or({
                if canonicalize_virtual_antiparticle {
                    particle.antiparticle
                } else {
                    edge.data.particle
                }
            });
            if (left, right) != edge.vertices {
                reversed_edges.insert(EdgeId(edge_id));
            }
            normalized_graph
                .add_edge(left, right, edge.directed, data)
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

    fn vertex_half_edge(
        model: &Model,
        particle: &Particle,
    ) -> Result<HalfEdge<EdgeColor>, GenerationError> {
        let base = if particle.is_antiparticle() {
            particle.antiparticle
        } else {
            model.particle_id(&particle.name)?
        };
        Ok(HalfEdge {
            direction: if model.particle_is_self_conjugate(model.particle_id(&particle.name)?) {
                None
            } else {
                Some(particle.is_antiparticle())
            },
            data: EdgeColor {
                particle: base,
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
        let he_graph: HedgeGraph<EdgeColor, ColoredNode, ()> = graph.clone().into();
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
            if let ColoredNode::External(external) = node {
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
            .iter()
            .enumerate()
            .filter_map(|(edge, data)| {
                let (source, target) = data.vertices;
                (matches!(graph.nodes()[source].data, ColoredNode::External(_))
                    || matches!(graph.nodes()[target].data, ColoredNode::External(_)))
                .then_some(EdgeId(edge))
            })
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
                let ColoredNode::Interaction(rule_id) = node else {
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
                    let particle = he_graph[[&hedge]].particle;
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
        let connected_components_before = graph
            .tadpoles(&externals)
            .expect("external nodes are distinct graph terminals with incident edges")
            .len()
            + 1;
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
    use crate::{CancellationToken, GraphGroupingOptions, NumeratorGrouping};
    use feynkit_model::Model;
    use symbolica::parser::ParseSettings;

    use super::*;
    use idenso::tensor::composition::TensorCompositionError;
    use idenso::tensor::inference::TensorInferenceError;
    use serde_json::json;
    use spenso::tensor_symbol;
    use std::array;

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
                    {"name":"f_prop","particle":"f","numerator":"PSlash(2,1)","denominator":"P^2"},
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
                    "color_structures":["1","1"],
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

    fn edge(model: &Model, pdg: i64) -> EdgeColor {
        EdgeColor {
            particle: particle(model, pdg),
            direction: None,
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
                .map(FeynmanDiagram::id)
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
                assert_eq!(
                    function_count(
                        &NumeratorInterface::export(diagram.projector()).0,
                        wavefunction
                    ),
                    1
                );
            }
            assert_eq!(
                function_count(
                    &NumeratorInterface::export(diagram.projector()).0,
                    feynkit_graph::symbols::hedge_index()
                ),
                4
            );
        }

        let supplied = Atom::mul_many(
            [
                (1, ExternalState::Incoming),
                (-1, ExternalState::Incoming),
                (1, ExternalState::Outgoing),
                (-1, ExternalState::Outgoing),
            ]
            .into_iter()
            .enumerate()
            .map(|(index, (pdg, state))| {
                generator
                    .external_wavefunction(state, particle(&generator.model, pdg), index, index)
                    .unwrap()
            }),
        );
        let overridden = generator
            .generate(
                &process,
                &options.clone().projector(supplied),
                GenerationType::Amplitude,
            )
            .unwrap();
        assert_eq!(overridden.diagrams.len(), generated.diagrams.len());
        for (automatic, supplied) in generated.diagrams.iter().zip(&overridden.diagrams) {
            assert_eq!(automatic.id(), supplied.id());
            assert_eq!(automatic.projector(), supplied.projector());
            assert_eq!(automatic.numerator(), supplied.numerator());
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
                .all(|diagram| diagram.projector().expression() == &Atom::one())
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
        let old_graph: HedgeGraph<EdgeColor, ColoredNode, ()> = graph.clone().into();
        let mut old_source = Vec::new();
        let mut old_target = Vec::new();
        for (node, _, color) in old_graph.iter_nodes() {
            let ColoredNode::External(external) = color else {
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
            .cut_partitions(&graph, &model, &GenerationOptions::default())
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
        let outgoing_0 = builder.add_generation_external("out0", 0, ExternalState::Outgoing, 0);
        let outgoing_1 = builder.add_generation_external("out1", 1, ExternalState::Outgoing, 1);
        let vertices = (0..8)
            .map(|vertex| {
                builder.add_vertex(DiagramVertex {
                    name: format!("v{vertex}"),
                    interaction: None,

                    numerator: SymbolicTensor::infer(Atom::one()).unwrap(),
                })
            })
            .collect::<Vec<_>>();
        assert_eq!(
            builder
                .add_edge(outgoing_0, vertices[7], DiagramEdge::new(particle, false),)
                .unwrap(),
            EdgeId(0),
        );
        assert_eq!(
            builder
                .add_edge(outgoing_1, vertices[7], DiagramEdge::new(particle, false),)
                .unwrap(),
            EdgeId(1),
        );
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

    #[test]
    fn supplied_vector_projector_uses_external_labels_before_topology_finalization() {
        let generator = Generator::new(fermion_model());
        for reversed in [false, true] {
            let mut graph = Graph::new();
            let left = graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V")));
            let right = graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V")));
            graph
                .add_edge(left, right, false, edge(&generator.model, 22))
                .unwrap();
            let mut supplied = Atom::one();
            for (index, state, attached) in [
                (7, ExternalState::Incoming, left),
                (3, ExternalState::Outgoing, right),
            ] {
                let external = graph.add_node(ColoredNode::External(ExternalNode {
                    index,
                    state,
                    particle: particle(&generator.model, 22),
                }));
                let (source, target) = if reversed {
                    (attached, external)
                } else {
                    (external, attached)
                };
                graph
                    .add_edge(source, target, false, edge(&generator.model, 22))
                    .unwrap();
                supplied *= generator
                    .external_wavefunction(state, particle(&generator.model, 22), index, index)
                    .unwrap();
            }
            assert_eq!(
                generator
                    .external_state_projector(&graph, Some(&supplied))
                    .unwrap(),
                generator.external_state_projector(&graph, None).unwrap(),
            );
            let invalid = feynkit_graph::symbols::epsilon().call(999usize);
            assert!(matches!(
                generator.external_state_projector(&graph, Some(&invalid)),
                Err(GenerationError::UnknownProjectorExternalIndex { index: 999 })
            ));
        }
    }

    #[test]
    fn finalizes_vector_projectors_and_leaves_scalar_states_trivial() {
        let generator = Generator::new(fermion_model());
        let mut vector_graph = Graph::new();
        let incoming = vector_graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Incoming,
            particle: particle(&generator.model, 22),
        }));
        let left = vector_graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V")));
        let right = vector_graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V")));
        let outgoing = vector_graph.add_node(ColoredNode::External(ExternalNode {
            index: 1,
            state: ExternalState::Outgoing,
            particle: particle(&generator.model, 22),
        }));
        vector_graph
            .add_edge(incoming, left, false, edge(&generator.model, 22))
            .unwrap();
        vector_graph
            .add_edge(right, outgoing, false, edge(&generator.model, 22))
            .unwrap();
        let vector_projector = generator
            .external_state_projector(&vector_graph, None)
            .unwrap();
        let vector_projector = NumeratorInterface::export(&vector_projector).0;
        assert_eq!(
            function_count(&vector_projector, feynkit_graph::symbols::epsilon()),
            1
        );
        assert_eq!(
            function_count(&vector_projector, feynkit_graph::symbols::epsilonbar()),
            1
        );
        assert_eq!(
            function_count(&vector_projector, feynkit_graph::symbols::hedge_index()),
            2
        );

        let mut scalar_graph = Graph::new();
        let scalar = scalar_graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Incoming,
            particle: particle(&generator.model, 25),
        }));
        let internal =
            scalar_graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V")));
        scalar_graph
            .add_edge(scalar, internal, false, edge(&generator.model, 25))
            .unwrap();
        assert_eq!(
            generator
                .external_state_projector(&scalar_graph, None)
                .unwrap()
                .expression(),
            &Atom::one()
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
                assert_eq!(serial.source_id, parallel.source_id);
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

    #[test]
    fn single_thread_generation_stays_on_the_caller() {
        let caller = std::thread::current().id();
        let observed = Arc::new(Mutex::new(Vec::new()));
        let sink = observed.clone();
        let options = GenerationOptions::default()
            .with_loop_count(1, 1)
            .unwrap()
            .threads(1)
            .filter_zero_color(true)
            .progress(move |progress| {
                if progress.completed > 0 {
                    sink.lock()
                        .unwrap()
                        .push((progress.stage.to_owned(), std::thread::current().id()));
                }
                crate::GenerationControl::Continue
            });
        let generator = Generator::new(Model::yang_mills());
        let process = Process::new(["g"], ["g"]);
        let first = generator
            .generate(&process, &options, GenerationType::Amplitude)
            .unwrap();
        let second = generator
            .generate(&process, &options, GenerationType::Amplitude)
            .unwrap();
        assert!(first.report.completed && second.report.completed);
        first.validate_groups().unwrap();
        second.validate_groups().unwrap();
        assert!(!first.diagrams.is_empty());
        assert_eq!(first.diagrams.len(), second.diagrams.len());
        for (left, right) in first.diagrams.iter().zip(&second.diagrams) {
            assert_eq!(left.to_dot().unwrap(), right.to_dot().unwrap());
        }
        let observed = observed.lock().unwrap();
        for stage in ["interactions", "numerators", "filter_zero_color"] {
            assert!(
                observed.iter().any(|(actual, _)| actual == stage),
                "missing {stage}"
            );
        }
        assert!(observed.iter().all(|(_, thread)| *thread == caller));
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn color_filter_uses_the_requested_generation_pool() {
        for threads in [2, 4] {
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
    fn directed_and_undirected_half_edges_use_the_interaction_key_order() {
        let model = fermion_model();
        let options =
            GenerationOptions::default().with_graph_filter(GenerationFilter::VertexAllow(vec![
                VertexSelector::Name("V".to_owned()),
            ]));
        let options = options.resolve_selectors(&model).unwrap();
        let signatures = InteractionSignatures::new(&model, &options).unwrap();
        let mut signature = vec![
            EdgeColor {
                particle: particle(&model, 1),
                direction: Some(false),
            },
            EdgeColor {
                particle: particle(&model, 1),
                direction: Some(true),
            },
            EdgeColor {
                particle: particle(&model, 22),
                direction: None,
            },
        ];
        signature.sort();

        assert_eq!(
            signatures.rules(&signature).unwrap(),
            &[vertex(&model, "V")]
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
        let fermion = EdgeColor {
            particle: particle(&model, 1),
            direction: None,
        };
        graph
            .add_edge(particle_vertex, left, true, fermion)
            .unwrap();
        graph
            .add_edge(
                right,
                left,
                true,
                EdgeColor {
                    particle: particle(&model, -1),
                    direction: None,
                },
            )
            .unwrap();
        graph.add_edge(right, antiparticle, true, fermion).unwrap();

        let normalized = process.normalize_fermion_flows(&graph, &model).unwrap();
        assert_eq!(
            normalized.graph.edges()[0].vertices,
            (particle_vertex, left)
        );
        assert_eq!(normalized.graph.edges()[1].vertices, (left, right));
        assert_eq!(normalized.graph.edges()[2].vertices, (right, antiparticle));
        assert_eq!(
            normalized
                .graph
                .edges()
                .iter()
                .map(|edge| edge.data.particle)
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
        .normalize_fermion_flows(&graph, &model)
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
        let flow = EdgeColor {
            particle: particle(&model, 1),
            direction: None,
        };
        graph.add_edge(ghost, left, true, flow).unwrap();
        graph.add_edge(right, left, true, flow).unwrap();
        graph.add_edge(right, antighost, true, flow).unwrap();

        let normalized = process.normalize_fermion_flows(&graph, &model).unwrap();
        assert_eq!(normalized.graph.edges()[0].vertices, (ghost, left));
        assert_eq!(normalized.graph.edges()[1].vertices, (left, right));
        assert_eq!(normalized.graph.edges()[2].vertices, (right, antighost));
        // The all-incoming order [ghost, antighost] is one transposition
        // from the canonical [antighost, ghost] order, just as for spinors.
        assert!(normalized.signs.external_ordering_negative);
        assert!(!normalized.signs.antifermion_spin_sum_negative);
    }

    #[test]
    fn cp_symmetrization_groups_mirrored_real_scalar_forward_graphs() {
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
                "name":"V_4_SCALAR_0000", "particles":["scalar_0","scalar_0","scalar_0","scalar_0"],
                "color_structures":["1"], "lorentz_structures":["SCALAR_4_LORENTZ_STRUCTURE"],
                "couplings":[["SCALAR_COUPLING"]]
            }));
        let model: Model = serde_json::from_value(definition).unwrap();
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
                        .canonicalize_numerator_graph(graph, &model, &NumeratorGrouping::None)
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
        let charged = EdgeColor {
            particle: particle(&model, 1),
            direction: None,
        };
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
        for node in graph.nodes() {
            mirror.add_node(node.data.clone());
        }
        for edge in graph.edges() {
            let mut data = edge.data;
            data.particle = model.particle_by_id(data.particle).unwrap().antiparticle;
            mirror
                .add_edge(
                    external_swap
                        .get(&edge.vertices.0)
                        .copied()
                        .unwrap_or(edge.vertices.0),
                    external_swap
                        .get(&edge.vertices.1)
                        .copied()
                        .unwrap_or(edge.vertices.1),
                    edge.directed,
                    data,
                )
                .unwrap();
        }

        let grouping = NumeratorGrouping::UpToScalar(GraphGroupingOptions::default());
        let canonical = process
            .canonicalize_numerator_graph(&graph, &model, &grouping)
            .unwrap();
        let canonical_mirror = process
            .canonicalize_numerator_graph(&mirror, &model, &grouping)
            .unwrap();
        assert_eq!(canonical, canonical_mirror);

        let internal_edge = canonical
            .edges()
            .iter()
            .enumerate()
            .find(|(_, edge)| {
                matches!(
                    canonical.nodes()[edge.vertices.0].data,
                    ColoredNode::Interaction(_)
                ) && matches!(
                    canonical.nodes()[edge.vertices.1].data,
                    ColoredNode::Interaction(_)
                )
            })
            .map(|(edge, _)| edge)
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
                    .particle_by_id(canonical.edges()[internal_edge].data.particle)
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
        let localized = NumeratorInstantiation {
            owner: NumeratorOwner::Edge(internal_edge),
            legs: &legs,
        }
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
                .canonicalize_numerator_graph(&graph, &model, &grouping)
                .unwrap();
            let pdgs = canonical
                .edges()
                .iter()
                .map(|edge| model.particle_by_id(edge.data.particle).unwrap().pdg_code)
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

        let repaired = ResolvedProcess::repair_cp_interaction_rules(graph, &model).unwrap();
        assert_eq!(
            repaired.nodes()[interaction].data,
            ColoredNode::Interaction(vertex(&model, "V_117"))
        );
        let rule = model.vertex_rule("V_117").unwrap();
        assert_eq!(
            Generator::new(model.clone())
                .interaction_half_edges(&repaired, interaction, rule)
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
        let repaired = ResolvedProcess::repair_cp_interaction_rules(graph, &model).unwrap();
        assert_eq!(
            repaired.nodes()[interaction].data,
            ColoredNode::Interaction(vertex(&model, "V_117"))
        );

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
            !generator
                .passes_interaction_filters(&repaired, 0, &allow_original)
                .unwrap(),
            "a CP repair must not bypass an allow-list"
        );
        assert!(
            !generator
                .passes_interaction_filters(&repaired, 0, &veto_repaired)
                .unwrap(),
            "a CP repair must not bypass a veto-list"
        );
        assert!(
            generator
                .passes_interaction_filters(&repaired, 0, &allow_repaired)
                .unwrap(),
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
            .add_edge(
                interaction,
                interaction,
                true,
                EdgeColor {
                    particle: particle(&model, -1),
                    direction: None,
                },
            )
            .unwrap();
        graph
            .add_edge(external, interaction, false, edge(&model, 22))
            .unwrap();

        let diagram = generator
            .to_diagram(
                "antiparticle-propagator".to_owned(),
                GenerationType::Amplitude,
                ColoredTopology {
                    graph,
                    symmetry: 1,
                    multiplicity: 1,
                    cut_partitions: Vec::new(),
                    topology_threshold_candidates: Vec::new(),
                },
                FermionSigns {
                    internal_loop_negative: true,
                    ..FermionSigns::default()
                },
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
            NumeratorInstantiation {
                owner: NumeratorOwner::Edge(0),
                legs: &legs,
            }
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

        let actual = NumeratorInterface::export(actual).0;
        assert_eq!(actual, exact);
        assert_ne!(actual, blanket_reversal);
        assert_ne!(actual, particle_fallback);
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
    fn vertex_orders_take_each_order_maximum_and_graph_orders_sum_vertices() {
        let model = fermion_model();
        let rule_orders = model.vertex_rule("V").unwrap().coupling_orders(&model);
        assert_eq!(
            rule_orders,
            BTreeMap::from([("QCD".to_owned(), 2), ("QED".to_owned(), 3)])
        );

        let generator = Generator::new(model);
        let mut graph = Graph::new();
        graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V")));
        graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V")));
        assert_eq!(
            generator.coupling_orders(&graph).unwrap(),
            BTreeMap::from([("QCD".to_owned(), 4), ("QED".to_owned(), 6)])
        );
    }

    #[test]
    fn interaction_numerator_keeps_named_couplings_and_lowers_tensors() {
        let mut definition = serde_json::to_value(fermion_model()).unwrap();
        for particle in definition["particles"].as_array_mut().unwrap() {
            particle["color"] = serde_json::json!(match particle["name"].as_str() {
                Some("f") => 3,
                Some("f~") => -3,
                Some("a") => 8,
                _ => 1,
            });
        }
        // Both couplings multiply the same open color interface. A scalar
        // color term cannot be added to a fundamental-adjoint tensor term.
        definition["vertex_rules"][0]["color_structures"] =
            serde_json::json!(["T(3,2,1)", "T(3,2,1)"]);
        let model = Model::from_json(&definition.to_string()).unwrap();
        let generator = Generator::new(model.clone());
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(7),
            legs: &legs[..3],
        };
        let numerator = generator
            .interaction_numerator(model.vertex_rule_id("V").unwrap(), &instantiation)
            .unwrap();

        let numerator = NumeratorInterface::export(&numerator).0.to_plain_string();
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
    fn interaction_numerator_rejects_incompatible_color_summands() {
        let mut definition = serde_json::to_value(fermion_model()).unwrap();
        definition["vertex_rules"][0]["color_structures"] = json!(["1", "T(3,2,1)"]);
        let model = Model::from_json(&definition.to_string()).unwrap();
        let generator = Generator::new(model.clone());
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(7),
            legs: &legs[..3],
        };
        let error = generator
            .interaction_numerator(model.vertex_rule_id("V").unwrap(), &instantiation)
            .unwrap_err();

        assert!(matches!(
            &error,
            GenerationError::Diagram(DiagramError::Tensor(_))
        ));
        assert!(
            error
                .to_string()
                .contains("tensor summands do not have compatible tensor interfaces")
        );
        assert!(generator.numerators.interactions.read().unwrap().is_empty());
        assert!(
            generator
                .numerators
                .interaction_templates
                .read()
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn interaction_templates_reuse_admission_and_localize_shifted_bound_indices() {
        let mut definition = serde_json::to_value(scalar_model()).unwrap();
        definition["lorentz_structures"][0]["structure"] =
            json!("Metric(idx(2,1),dummy(1))*P(dummy(1),2)");
        let model = Model::from_json(&definition.to_string()).unwrap();
        let rule = model.vertex_rule_id("V1").unwrap();
        let generator = Generator::new(model);
        let mut legs = tensor_test_legs()[..3].to_vec();
        let mut previous = None;
        for owner in [7, 12] {
            for leg in &mut legs {
                leg.edge += 13;
                leg.hedge += 19;
                leg.color = 1;
            }
            let instantiation = NumeratorInstantiation {
                owner: NumeratorOwner::Vertex(owner),
                legs: &legs,
            };
            let actual = generator
                .interaction_numerator(rule, &instantiation)
                .unwrap();
            let expected = NumeratorInterface::admit(
                generator
                    .interaction_expression(rule, &instantiation)
                    .unwrap(),
                None,
            )
            .unwrap();
            assert!(actual.same_value(&expected));
            assert_eq!(actual.rank(), 1);
            assert_eq!(
                actual.structure().logical_slots()[0].aind,
                PartialIndex::Explicit(AbstractIndex::Named(
                    symbols::hedge_index().into(),
                    legs[0].hedge,
                    2,
                ))
            );
            let templates = generator.numerators.interaction_templates.read().unwrap();
            assert_eq!(templates.len(), 1);
            let template = templates.values().next().unwrap().clone();
            if let Some(previous) = previous {
                assert!(Arc::ptr_eq(&previous, &template));
            }
            previous = Some(template);
        }
        assert_eq!(generator.numerators.interactions.read().unwrap().len(), 2);
    }

    #[test]
    fn interaction_template_cache_retains_dummy_metadata_across_owners() {
        let mut definition = serde_json::to_value(scalar_model()).unwrap();
        definition["lorentz_structures"][0]["structure"] = json!("P(dummy(1),1)*P(dummy(1),2)");
        let model = Model::from_json(&definition.to_string()).unwrap();
        let rule = model.vertex_rule_id("V1").unwrap();
        let generator = Generator::new(model);
        let mut legs = tensor_test_legs()[..3].to_vec();
        for leg in &mut legs {
            leg.color = 1;
        }
        let head = tensor_symbol!("generation_template_cache_probe");
        let mut numerators = Vec::new();
        for owner in [7, 12] {
            let instantiation = NumeratorInstantiation {
                owner: NumeratorOwner::Vertex(owner),
                legs: &legs,
            };
            let numerator = generator
                .interaction_numerator(rule, &instantiation)
                .unwrap();
            let index = AbstractIndex::Named(symbols::vertex_index().into(), owner, 1);
            let probe = SymbolicTensor::infer(
                FunctionBuilder::new(head)
                    .add_arg(Minkowski {}.new_rep(4).to_symbolic([index.to_atom()]))
                    .finish(),
            )
            .unwrap();
            assert!(matches!(
                numerator.multiply(&probe),
                Err(TensorCompositionError::InvalidExplicitMultiplicity {
                    index: collision,
                    occurrences: 3,
                }) if collision == index
            ));
            // Consuming a returned clone must not establish metadata which the
            // retained cache entry lacks and would have to rediscover later.
            let cached = generator
                .interaction_numerator(rule, &instantiation)
                .unwrap();
            assert!(numerator.same_reduction_state(&cached));
            assert!(matches!(
                cached.multiply(&probe),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            ));
            numerators.push(cached);
        }
        assert!(numerators[0].multiply(&numerators[1]).is_ok());
        assert!(numerators[0].multiply(&numerators[0]).is_err());
        assert_eq!(
            generator
                .numerators
                .interaction_templates
                .read()
                .unwrap()
                .len(),
            1
        );
    }

    #[test]
    fn cancelled_alias_templates_keep_hidden_dummy_labels_and_cached_metadata() {
        let mut definition = serde_json::to_value(scalar_model()).unwrap();
        definition["lorentz_structures"][0]["structure"] =
            json!("Metric(idx(2,3),dummy(1))*(P(dummy(1),1)+P(dummy(1),2))");
        let model = Model::from_json(&definition.to_string()).unwrap();
        let rule = model.vertex_rule_id("V1").unwrap();
        let generator = Generator::new(model);
        let mut legs = tensor_test_legs()[..3].to_vec();
        for leg in &mut legs {
            leg.color = 1;
        }
        legs[1].edge = legs[0].edge;
        let head = tensor_symbol!("generation_cancelled_template_probe");
        for owner in [4, 9] {
            let instantiation = NumeratorInstantiation {
                owner: NumeratorOwner::Vertex(owner),
                legs: &legs,
            };
            let zero = generator
                .interaction_numerator(rule, &instantiation)
                .unwrap();
            let cached = generator
                .interaction_numerator(rule, &instantiation)
                .unwrap();
            assert!(zero.is_zero());
            assert_eq!(zero.rank(), 1);
            assert!(zero.same_reduction_state(&cached));
            let index = AbstractIndex::Named(symbols::vertex_index().into(), owner, 1);
            let probe = |index: AbstractIndex| {
                SymbolicTensor::infer(
                    FunctionBuilder::new(head)
                        .add_arg(Minkowski {}.new_rep(4).to_symbolic([index.to_atom()]))
                        .finish(),
                )
                .unwrap()
            };
            for numerator in [zero, cached] {
                assert!(matches!(
                    numerator.multiply(&probe(index)),
                    Err(TensorCompositionError::InvalidExplicitMultiplicity {
                        index: collision,
                        occurrences: 3,
                    }) if collision == index
                ));
                let other = AbstractIndex::Named(symbols::vertex_index().into(), owner + 1, 1);
                let product = numerator.multiply(&probe(other)).unwrap();
                assert!(product.is_zero());
                assert_eq!(product.rank(), 2);
            }
        }
        let templates = generator.numerators.interaction_templates.read().unwrap();
        assert_eq!(templates.len(), 2);
        let dummy = AbstractIndex::Named(symbols::vertex_index().into(), 0, 1);
        for template in templates.values() {
            let NumeratorTemplate::Admitted { labels, .. } = &**template else {
                panic!("explicit rule must retain its admitted template");
            };
            assert!(labels.iter().any(|(index, _)| *index == dummy));
        }
    }

    #[test]
    fn interaction_templates_preserve_rank_when_incident_momenta_cancel() {
        let mut definition = serde_json::to_value(scalar_model()).unwrap();
        definition["lorentz_structures"][0]["structure"] = json!("P(idx(2,1),1)+P(idx(2,1),2)");
        let model = Model::from_json(&definition.to_string()).unwrap();
        let rule = model.vertex_rule_id("V1").unwrap();
        let generator = Generator::new(model);
        let mut legs = tensor_test_legs()[..3].to_vec();
        for leg in &mut legs {
            leg.color = 1;
        }
        legs[1].edge = legs[0].edge;
        for owner in [4, 9] {
            let instantiation = NumeratorInstantiation {
                owner: NumeratorOwner::Vertex(owner),
                legs: &legs,
            };
            let zero = generator
                .interaction_numerator(rule, &instantiation)
                .unwrap();
            assert!(zero.is_zero());
            assert_eq!(zero.rank(), 1);
            assert_eq!(
                zero.structure().logical_slots()[0].aind,
                PartialIndex::Explicit(AbstractIndex::Named(
                    symbols::hedge_index().into(),
                    legs[0].hedge,
                    2,
                ))
            );
            legs[0].edge += 6;
            legs[1].edge = legs[0].edge;
            for leg in &mut legs {
                leg.hedge += 7;
            }
        }
        assert_eq!(
            generator
                .numerators
                .interaction_templates
                .read()
                .unwrap()
                .len(),
            2
        );
    }

    #[test]
    fn propagator_templates_reuse_reversed_spin_and_color_interfaces() {
        let mut definition = serde_json::to_value(fermion_model()).unwrap();
        for particle in definition["particles"].as_array_mut().unwrap() {
            particle["color"] = json!(match particle["name"].as_str() {
                Some("f") => 3,
                Some("f~") => -3,
                _ => 1,
            });
        }
        let model = Model::from_json(&definition.to_string()).unwrap();
        let generator = Generator::new(model.clone());
        for pdg in [1, -1] {
            let particle_id = particle(&model, pdg);
            let particle = model.particle_by_id(particle_id).unwrap();
            let (propagator, reversed) = generator.propagator_numerator(particle).unwrap();
            for (edge, hedges) in [(5, [11, 17]), (8, [27, 21])] {
                let actual = generator.edge_numerator(particle_id, edge, hedges).unwrap();
                let flows = if reversed {
                    [Flow::Sink, Flow::Source]
                } else {
                    [Flow::Source, Flow::Sink]
                };
                let legs = flows.map(|flow| NumeratorHalfEdge {
                    edge,
                    hedge: hedges[usize::from(flow == Flow::Sink)],
                    flow,
                    spin: particle.spin,
                    color: particle.color,
                });
                let instantiation = NumeratorInstantiation {
                    owner: NumeratorOwner::Edge(edge),
                    legs: &legs,
                };
                let expected = NumeratorInterface::admit(
                    instantiation
                        .propagator_expression(&propagator, particle.color)
                        .unwrap(),
                    None,
                )
                .unwrap();
                assert_eq!(actual.expression(), expected.expression());
                assert_eq!(
                    actual.structure().logical_slots(),
                    expected.structure().logical_slots(),
                );
            }
        }
        assert_eq!(
            generator
                .numerators
                .propagator_templates
                .read()
                .unwrap()
                .len(),
            2
        );
    }

    #[test]
    fn interaction_template_order_matches_actual_owner_admission() {
        for expression in [
            "Metric(1,2)",
            "P(idx(1,1),1)*P(idx(1,2),2)",
            "Metric(1,2)*(P(3,1)-P(3,2))",
        ] {
            let mut definition = serde_json::to_value(scalar_model()).unwrap();
            definition["lorentz_structures"][0]["structure"] = json!(expression);
            let model = Model::from_json(&definition.to_string()).unwrap();
            let rule = model.vertex_rule_id("V1").unwrap();
            let generator = Generator::new(model);
            for (owner, permutation) in [[0, 1, 2], [1, 0, 2], [2, 0, 1], [2, 1, 0]]
                .into_iter()
                .enumerate()
            {
                let legs = permutation.map(|index| NumeratorHalfEdge {
                    edge: [17, 3, 23][index],
                    hedge: [41, 11, 29][index],
                    flow: Flow::Sink,
                    spin: 1,
                    color: 1,
                });
                let instantiation = NumeratorInstantiation {
                    owner: NumeratorOwner::Vertex(owner),
                    legs: &legs,
                };
                let actual = generator
                    .interaction_numerator(rule, &instantiation)
                    .unwrap();
                let expected = NumeratorInterface::admit(
                    generator
                        .interaction_expression(rule, &instantiation)
                        .unwrap(),
                    None,
                )
                .unwrap();
                assert_eq!(
                    actual.expression().as_view().get_data(),
                    expected.expression().as_view().get_data(),
                    "{expression}: {permutation:?}"
                );
                assert_eq!(
                    actual.structure().logical_slots(),
                    expected.structure().logical_slots(),
                    "{expression}: {permutation:?}"
                );
            }
            let templates = generator.numerators.interaction_templates.read().unwrap();
            assert_eq!(templates.len(), 1);
            assert!(matches!(
                &**templates.values().next().unwrap(),
                NumeratorTemplate::Admitted { .. }
            ));
        }
    }

    #[test]
    fn propagator_template_order_matches_actual_owner_admission() {
        for expression in ["Metric(1,2)", "P(idx(1,1))*P(idx(1,2))"] {
            let mut definition = serde_json::to_value(scalar_model()).unwrap();
            definition["propagators"][0]["numerator"] = json!(expression);
            let model = Model::from_json(&definition.to_string()).unwrap();
            let particle_id = particle(&model, 25);
            let generator = Generator::new(model);
            let particle = generator.model.particle_by_id(particle_id).unwrap();
            let (propagator, reversed) = generator.propagator_numerator(particle).unwrap();
            assert!(!reversed);
            for (edge, hedges) in [(5, [11, 17]), (8, [27, 21])] {
                let actual = generator.edge_numerator(particle_id, edge, hedges).unwrap();
                let legs = [Flow::Source, Flow::Sink].map(|flow| NumeratorHalfEdge {
                    edge,
                    hedge: hedges[usize::from(flow == Flow::Sink)],
                    flow,
                    spin: particle.spin,
                    color: particle.color,
                });
                let instantiation = NumeratorInstantiation {
                    owner: NumeratorOwner::Edge(edge),
                    legs: &legs,
                };
                let expected = NumeratorInterface::admit(
                    instantiation
                        .propagator_expression(&propagator, particle.color)
                        .unwrap(),
                    None,
                )
                .unwrap();
                assert_eq!(
                    actual.expression().as_view().get_data(),
                    expected.expression().as_view().get_data(),
                    "{expression}: {hedges:?}"
                );
                assert_eq!(
                    actual.structure().logical_slots(),
                    expected.structure().logical_slots(),
                    "{expression}: {hedges:?}"
                );
            }
            let templates = generator.numerators.propagator_templates.read().unwrap();
            assert_eq!(templates.len(), 1);
            assert!(matches!(
                &**templates.values().next().unwrap(),
                NumeratorTemplate::Admitted { .. }
            ));
        }
    }

    #[test]
    fn interaction_scoped_pairs_do_not_reorder_public_ports() {
        let _ = tensor_symbol!("UFO::generation_scoped_order_tensor");
        let _ = Minkowski {}.new_rep(4).to_symbolic([]);
        for scalar in [
            "(P(1,1)*P(1,2)+P(1,2)*P(1,3))",
            "spenso::bracket(P(1,1),P(1,2))^2",
        ] {
            let expression = format!(
                "{scalar}*generation_scoped_order_tensor(spenso::mink(4,idx(1,2)),spenso::mink(4,idx(1,1)))"
            );
            let mut definition = serde_json::to_value(scalar_model()).unwrap();
            definition["lorentz_structures"][0]["structure"] = json!(expression);
            let model = Model::from_json(&definition.to_string()).unwrap();
            let rule = model.vertex_rule_id("V1").unwrap();
            let generator = Generator::new(model);
            for hedges in [[0, 1, 2], [27, 11, 19]] {
                let legs = array::from_fn::<_, 3, _>(|index| NumeratorHalfEdge {
                    edge: index + 3,
                    hedge: hedges[index],
                    flow: Flow::Sink,
                    spin: 1,
                    color: 1,
                });
                let instantiation = NumeratorInstantiation {
                    owner: NumeratorOwner::Vertex(7),
                    legs: &legs,
                };
                let expected = NumeratorInterface::admit(
                    generator
                        .interaction_expression(rule, &instantiation)
                        .unwrap(),
                    None,
                )
                .unwrap();
                let actual = generator
                    .interaction_numerator(rule, &instantiation)
                    .unwrap();
                let expected_indices = expected
                    .structure()
                    .logical_slots()
                    .into_iter()
                    .map(|slot| slot.aind)
                    .collect::<Vec<_>>();
                assert_eq!(
                    expected_indices,
                    [hedges[1], hedges[0]].map(|hedge| {
                        PartialIndex::Explicit(AbstractIndex::Named(
                            symbols::hedge_index().into(),
                            hedge,
                            1,
                        ))
                    })
                );
                assert_eq!(
                    actual.expression().as_view().get_data(),
                    expected.expression().as_view().get_data(),
                    "{scalar}"
                );
                assert_eq!(
                    actual.structure().logical_slots(),
                    expected.structure().logical_slots(),
                    "{scalar}"
                );
            }
            let templates = generator.numerators.interaction_templates.read().unwrap();
            assert_eq!(templates.len(), 1);
            assert!(matches!(
                &**templates.values().next().unwrap(),
                NumeratorTemplate::Dynamic
            ));
        }
    }

    #[test]
    fn propagator_scoped_pairs_do_not_reorder_public_ports() {
        let _ = tensor_symbol!("UFO::generation_scoped_prop_tensor");
        let _ = Minkowski {}.new_rep(4).to_symbolic([]);
        for scalar in ["(P(1,1)*P(1,2)+1)", "spenso::bracket(P(1,1),P(1,2))^2"] {
            let expression = format!(
                "{scalar}*generation_scoped_prop_tensor(spenso::mink(4,idx(1,2)),spenso::mink(4,idx(1,1)))"
            );
            let mut definition = serde_json::to_value(scalar_model()).unwrap();
            definition["propagators"][0]["numerator"] = json!(expression);
            let model = Model::from_json(&definition.to_string()).unwrap();
            let particle_id = particle(&model, 25);
            let generator = Generator::new(model);
            let particle = generator.model.particle_by_id(particle_id).unwrap();
            let (propagator, reversed) = generator.propagator_numerator(particle).unwrap();
            assert!(!reversed);
            for (edge, hedges) in [(3, [0, 1]), (5, [27, 11])] {
                let legs = [Flow::Source, Flow::Sink].map(|flow| NumeratorHalfEdge {
                    edge,
                    hedge: hedges[usize::from(flow == Flow::Sink)],
                    flow,
                    spin: particle.spin,
                    color: particle.color,
                });
                let instantiation = NumeratorInstantiation {
                    owner: NumeratorOwner::Edge(edge),
                    legs: &legs,
                };
                let expected = NumeratorInterface::admit(
                    instantiation
                        .propagator_expression(&propagator, particle.color)
                        .unwrap(),
                    None,
                )
                .unwrap();
                let actual = generator.edge_numerator(particle_id, edge, hedges).unwrap();
                let expected_indices = expected
                    .structure()
                    .logical_slots()
                    .into_iter()
                    .map(|slot| slot.aind)
                    .collect::<Vec<_>>();
                assert_eq!(
                    expected_indices,
                    [hedges[1], hedges[0]].map(|hedge| {
                        PartialIndex::Explicit(AbstractIndex::Named(
                            symbols::hedge_index().into(),
                            hedge,
                            1,
                        ))
                    })
                );
                assert_eq!(
                    actual.expression().as_view().get_data(),
                    expected.expression().as_view().get_data(),
                    "{scalar}"
                );
                assert_eq!(
                    actual.structure().logical_slots(),
                    expected.structure().logical_slots(),
                    "{scalar}"
                );
            }
            let templates = generator.numerators.propagator_templates.read().unwrap();
            assert_eq!(templates.len(), 1);
            assert!(matches!(
                &**templates.values().next().unwrap(),
                NumeratorTemplate::Dynamic
            ));
        }
    }

    #[test]
    fn callback_rules_are_constructed_in_the_actual_owner_frame() {
        let callback = symbol!(
            "UFO::generation_owner_callback",
            norm = |node, out| {
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(index)) = function.iter().next()
                    && index.get_symbol() == symbols::hedge_index()
                {
                    **out = index.iter().next().expect("owner ID").to_owned();
                }
            }
        );
        let mut definition = serde_json::to_value(scalar_model()).unwrap();
        definition["lorentz_structures"][0]["structure"] =
            json!("generation_owner_callback(idx(1,1))");
        let model = Model::from_json(&definition.to_string()).unwrap();
        let rule = model.vertex_rule_id("V1").unwrap();
        let generator = Generator::new(model);
        let mut legs = tensor_test_legs()[..3].to_vec();
        for owner in [3, 6] {
            let instantiation = NumeratorInstantiation {
                owner: NumeratorOwner::Vertex(owner),
                legs: &legs,
            };
            let actual = generator
                .interaction_numerator(rule, &instantiation)
                .unwrap();
            assert_eq!(
                actual.expression(),
                &(Atom::num(legs[0].hedge) * Atom::var(symbol!("UFO::GC1")))
            );
            assert_eq!(function_count(actual.expression(), callback), 0);
            legs[0].hedge += 5;
        }
        assert!(
            generator
                .numerators
                .interaction_templates
                .read()
                .unwrap()
                .values()
                .all(|template| matches!(&**template, NumeratorTemplate::Dynamic))
        );
    }

    #[test]
    fn compact_rule_powers_allocate_independent_dummies_for_each_owner() {
        let _ = Minkowski {}.new_rep(4).to_symbolic([]);
        let _ = tensor_symbol!("UFO::generation_template_auto_tensor");
        for expression in [
            "generation_template_auto_tensor(spenso::mink(4))^2",
            "P(spenso::mink(4))^4",
        ] {
            let mut definition = serde_json::to_value(scalar_model()).unwrap();
            definition["lorentz_structures"][0]["structure"] = json!(expression);
            let model = Model::from_json(&definition.to_string()).unwrap();
            let rule = model.vertex_rule_id("V1").unwrap();
            let generator = Generator::new(model);
            let mut legs = tensor_test_legs()[..3].to_vec();
            let mut numerators = Vec::new();
            let mut owner_indices = Vec::new();
            for owner in [7, 12] {
                let instantiation = NumeratorInstantiation {
                    owner: NumeratorOwner::Vertex(owner),
                    legs: &legs,
                };
                let numerator = generator
                    .interaction_numerator(rule, &instantiation)
                    .unwrap();
                assert_eq!(numerator.rank(), 0);
                assert!(
                    numerator.same_value(
                        &generator
                            .interaction_numerator(rule, &instantiation)
                            .unwrap(),
                    )
                );
                let mut slots = SlotMatcher::default();
                let mut indices = HashSet::new();
                numerator.expression().visitor(&mut |node| {
                    if let SlotMatch::Explicit(slot) = slots.classify(node) {
                        indices.insert(AbstractIndex::try_from(slot.index()).unwrap());
                    }
                    true
                });
                assert!(
                    !indices.is_empty(),
                    "{expression} must materialize a dummy pair"
                );
                owner_indices.push(indices);
                numerators.push(numerator);
                for leg in &mut legs {
                    leg.edge += 13;
                    leg.hedge += 19;
                }
            }
            assert!(
                owner_indices[0].is_disjoint(&owner_indices[1]),
                "{expression}"
            );
            let product = numerators[0].multiply(&numerators[1]).unwrap();
            assert_eq!(product.rank(), 0);
            SymbolicTensor::validate_interface(product.expression(), product.structure()).unwrap();
            assert!(
                generator
                    .numerators
                    .interaction_templates
                    .read()
                    .unwrap()
                    .values()
                    .all(|template| matches!(&**template, NumeratorTemplate::Dynamic))
            );
        }
    }

    #[test]
    fn template_eligibility_rejects_scoped_and_encoded_open_indices() {
        let scope = symbol!("generation_template_open_scope");
        let open = AbstractIndex::Open { owner: 3, axis: 0 };
        let scoped = open.scoped(scope).to_atom();
        let encoded = CookSettings::reversible()
            .cook_function(scoped.as_fun_view().unwrap())
            .unwrap();
        let head = tensor_symbol!("generation_template_open_tensor");
        for index in [open.to_atom(), scoped, encoded] {
            let expression = FunctionBuilder::new(head)
                .add_arg(Minkowski {}.new_rep(4).to_symbolic([index]))
                .finish();
            assert!(!NumeratorTemplate::supports(&expression));
        }
    }

    #[test]
    fn bundled_rule_ports_preserve_actual_owner_rejection() {
        let head = tensor_symbol!("UFO::generation_template_bundle_tensor");
        let slot = Minkowski {}.new_rep(4).to_symbolic([test_atom("idx(1,1)")]);
        let expression = head.call((AIND_SYMBOLS.aind.call((slot,)),));
        let mut definition = serde_json::to_value(scalar_model()).unwrap();
        definition["lorentz_structures"][0]["structure"] = json!(expression.to_canonical_string());
        let model = Model::from_json(&definition.to_string()).unwrap();
        let rule = model.vertex_rule_id("V1").unwrap();
        let generator = Generator::new(model);
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(7),
            legs: &legs[..3],
        };
        let expected = NumeratorInterface::admit(
            generator
                .interaction_expression(rule, &instantiation)
                .unwrap(),
            None,
        )
        .unwrap_err();
        let actual = generator
            .interaction_numerator(rule, &instantiation)
            .unwrap_err();
        // Spenso's structural parser can read bundles, but typed admission
        // rejects the untagged aind head. Template selection must preserve
        // that actual-owner rejection rather than publish a trusted carrier.
        assert!(matches!(
            &expected,
            DiagramError::Tensor(TensorInferenceError::Invalid(message))
                if message.contains("carries structural ports but is not tagged as a tensor")
        ));
        let GenerationError::Diagram(actual) = actual else {
            panic!("expected the existing typed-admission rejection");
        };
        assert_eq!(actual.to_string(), expected.to_string());
        assert!(generator.numerators.interactions.read().unwrap().is_empty());
        let templates = generator.numerators.interaction_templates.read().unwrap();
        assert_eq!(templates.len(), 1);
        assert!(
            templates
                .values()
                .all(|template| matches!(&**template, NumeratorTemplate::Dynamic))
        );
    }

    #[test]
    fn opaque_scalar_metadata_with_nonindex_slot_payload_keeps_checked_construction() {
        let _ = symbol!("UFO::generation_opaque_metadata"; Scalar);
        let _ = Minkowski {}.new_rep(4).to_symbolic([]);
        let mut definition = serde_json::to_value(scalar_model()).unwrap();
        definition["lorentz_structures"][0]["structure"] =
            json!("generation_opaque_metadata(spenso::mink(4,x+y))");
        let model = Model::from_json(&definition.to_string()).unwrap();
        let rule = model.vertex_rule_id("V1").unwrap();
        let generator = Generator::new(model);
        let legs = tensor_test_legs();
        for owner in [3, 6] {
            let instantiation = NumeratorInstantiation {
                owner: NumeratorOwner::Vertex(owner),
                legs: &legs[..3],
            };
            let expected = NumeratorInterface::admit(
                generator
                    .interaction_expression(rule, &instantiation)
                    .unwrap(),
                None,
            )
            .unwrap();
            let actual = generator
                .interaction_numerator(rule, &instantiation)
                .unwrap();
            assert!(actual.same_value(&expected));
            assert_eq!(actual.rank(), 0);
        }
        assert!(
            generator
                .numerators
                .interaction_templates
                .read()
                .unwrap()
                .values()
                .all(|template| matches!(&**template, NumeratorTemplate::Dynamic))
        );
    }

    #[test]
    fn template_transport_retains_declared_order_scoped_labels_and_zero_ports() {
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let scope = symbol!("generation_template_scope");
        let labels = |hedges: [usize; 2]| {
            let plain = AbstractIndex::Named(symbols::hedge_index().into(), hedges[0], 2);
            let payload = symbol!("generation_template_label")
                .call((symbols::hedge_index().call((hedges[1], 3)),));
            let encoded = cooking
                .cook_function(payload.as_fun_view().unwrap())
                .unwrap();
            let scoped = AbstractIndex::try_from(encoded.as_view())
                .unwrap()
                .scoped(scope);
            [plain, scoped].map(|index| {
                Minkowski {}
                    .new_rep(4)
                    .cast::<LibraryRep>()
                    .slot::<PartialIndex<AbstractIndex>, _>(PartialIndex::Explicit(index))
            })
        };
        let original = labels([0, 1]);
        let head = tensor_symbol!("generation_template_order");
        let expression = FunctionBuilder::new(head)
            .add_args(original.map(|slot| match slot.aind {
                PartialIndex::Explicit(index) => {
                    slot.rep().slot::<AbstractIndex, _>(index).to_atom()
                }
                _ => unreachable!(),
            }))
            .finish();
        let source = SymbolicTensor::from_declared_parts(
            expression,
            PartialStructure::from_logical_slots([original[1], original[0]]),
        )
        .unwrap();
        let zero = source.with_checked_expression(Atom::Zero).unwrap();
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(7),
            legs: &legs[..2],
        };
        let target = labels([legs[0].hedge, legs[1].hedge]);
        for source in [source, zero] {
            let was_zero = source.is_zero();
            let template = NumeratorTemplate::admitted(source).unwrap();
            let mapped = template.instantiate(&instantiation).unwrap();
            assert_eq!(
                mapped.structure().logical_slots(),
                vec![target[1], target[0]]
            );
            assert_eq!(mapped.is_zero(), was_zero);
            SymbolicTensor::validate_interface(mapped.expression(), mapped.structure()).unwrap();
        }
    }

    #[test]
    fn lowers_ufo_lorentz_tensors_to_spenso_structures() {
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(7),
            legs: &legs,
        };
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
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Edge(30),
            legs: &legs[2..],
        };
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
    fn ufo_lorentz_epsilon_matches_aloha_and_native_gamma5() {
        use idenso::{CookMode, CookSettings, tensor::SymbolicTensor};
        use spenso::{
            network::{
                ExecutionResult, Sequential, SmallestDegree,
                library::symbolic::ExplicitKey,
                parsing::{ParseSettings as TensorParseSettings, StrictTensorFilter},
            },
            structure::abstract_index::{AIND_SYMBOLS, AbstractIndex},
            tensors::parametric::MixedTensor,
        };
        use spenso_hep_lib::{FUN_LIB, HepNet, hep_lib_atom};
        use symbolica::function;

        let mut legs = tensor_test_legs();
        for leg in &mut legs {
            leg.spin = 3;
        }
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(7),
            legs: &legs,
        };
        let lowered = instantiation
            .instantiate(&test_atom("Epsilon(1,2,3,4)"), NumeratorSector::Spin)
            .unwrap();
        let native_epsilon = test_atom(
            "spenso::epsilon(spenso::mink(4,gammalooprs::hedge(20,1)),spenso::mink(4,gammalooprs::hedge(41,1)),spenso::mink(4,gammalooprs::hedge(60,1)),spenso::mink(4,gammalooprs::hedge(81,1)))",
        );
        assert_eq!(lowered, -Atom::i() * &native_epsilon);
        assert_eq!(instantiation.lower_spin(lowered.clone()).unwrap(), lowered);
        assert_eq!(
            instantiation
                .instantiate(&test_atom("Epsilon(2,1,3,4)"), NumeratorSector::Spin)
                .unwrap(),
            -&lowered,
        );

        // Use the unchanged gamma5/gamma library as an exact component oracle,
        // not a second epsilon table. Idenso defines Tr(g5*g^a*g^b*g^c*g^d)
        // = 4*epsilon(a,b,c,d), with epsilon(0,1,2,3) = -i.
        let trace = instantiation
            .instantiate(
                &test_atom(
                    "Gamma5(-1,-2)*Gamma(1,-2,-3)*Gamma(2,-3,-4)*Gamma(3,-4,-5)*Gamma(4,-5,-1)",
                ),
                NumeratorSector::Spin,
            )
            .unwrap();
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let simplified = SymbolicTensor::infer(cooking.cook(trace.as_view()))
            .unwrap()
            .simplify_algebra(&idenso::tensor::AlgebraSettings {
                gamma: Some(idenso::dirac::GammaSimplifySettings::default()),
                ..Default::default()
            })
            .unwrap();
        assert_eq!(
            cooking.uncook(simplified.expression().as_view()),
            Atom::num(4) * &native_epsilon,
        );
        let component_expression = lowered
            .replace(native_epsilon.to_pattern())
            .with(&trace / Atom::num(4));
        let library = hep_lib_atom::<AbstractIndex, MixedTensor<f64, ExplicitKey<AbstractIndex>>>();
        let settings = TensorParseSettings::default()
            .with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
        // MadGraph/ALOHA L_Epsilon.give_parity gives (0,1,2,3) -> -1.
        // Include an odd permutation and a repeated index. Four Minkowski
        // selectors add a minus sign for the three spatial directions.
        for (coordinates, aloha_component) in [
            ([0, 1, 2, 3], -1),
            ([1, 0, 2, 3], 1),
            ([0, 2, 3, 1], -1),
            ([0, 1, 2, 2], 0),
        ] {
            let mut selected = component_expression.clone();
            for (coordinate, hedge) in coordinates.into_iter().zip([20, 41, 60, 81]) {
                selected *= function!(
                    ETS.delta,
                    function!(AIND_SYMBOLS.cind, coordinate),
                    test_atom(&format!("spenso::mink(4,gammalooprs::hedge({hedge},1))"))
                );
            }
            let selected = cooking.cook(selected.as_view());
            let mut network =
                HepNet::<AbstractIndex>::try_from_view(selected.as_view(), &library, &settings)
                    .unwrap();
            network
                .execute::<Sequential, SmallestDegree, _, _, _>(&library, &*FUN_LIB)
                .unwrap();
            let actual = match network.result_scalar().unwrap() {
                ExecutionResult::One => Atom::one(),
                ExecutionResult::Zero => Atom::Zero,
                ExecutionResult::Val(value) => value.into_owned(),
            };
            assert_eq!(actual, Atom::num(-aloha_component), "{coordinates:?}");
        }
    }

    #[test]
    fn lowers_ufo_epsilon_with_the_ufo_normalization() {
        use idenso::{
            CookMode, CookSettings,
            tensor::{AlgebraSettings, SymbolicTensor},
        };

        let legs: [NumeratorHalfEdge; 4] = std::array::from_fn(|i| NumeratorHalfEdge {
            edge: 10 * (i + 1),
            flow: Flow::Sink,
            hedge: 20 * (i + 1),
            spin: 3,
            color: 1,
        });
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(7),
            legs: &legs,
        };
        let lower = |structure: &str| {
            instantiation
                .instantiate(&test_atom(structure), NumeratorSector::Spin)
                .unwrap()
        };
        assert_eq!(
            lower("Epsilon(1,2,3,4)"),
            test_atom(
                "-1i*spenso::epsilon(spenso::mink(4,gammalooprs::hedge(20,1)),spenso::mink(4,gammalooprs::hedge(40,1)),spenso::mink(4,gammalooprs::hedge(60,1)),spenso::mink(4,gammalooprs::hedge(80,1)))"
            )
        );
        // ALOHA's explicit matrices satisfy
        // Tr[γ5 γ^μ γ^ν γ^ρ γ^σ] = 4i ε^{μνρσ} with its ε^{0123} = −1; the
        // lowered structures must agree after Idenso's trace simplification.
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let simplify = |structure: &str| {
            let tensor = SymbolicTensor::infer(cooking.cook(lower(structure).as_view())).unwrap();
            let simplified = tensor
                .simplify_algebra(&AlgebraSettings::hep())
                .unwrap()
                .into_expression();
            cooking.uncook(simplified.as_view())
        };
        let trace =
            simplify("Gamma5(-1,-2)*Gamma(1,-2,-3)*Gamma(2,-3,-4)*Gamma(3,-4,-5)*Gamma(4,-5,-1)");
        assert!(!trace.is_zero());
        assert_eq!(trace, simplify("4*1i*Epsilon(1,2,3,4)"));
        for malformed in ["Epsilon()", "Epsilon(1,2,3)", "Epsilon(1,2,3,4,-1)"] {
            assert!(matches!(
                instantiation.instantiate(&test_atom(malformed), NumeratorSector::Spin),
                Err(GenerationError::InvalidNumeratorTensor { .. })
            ));
        }
    }

    #[test]
    fn lowers_ufo_charge_conjugation_and_rejects_invalid_arity() {
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(7),
            legs: &legs,
        };
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
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(7),
            legs: &legs,
        };
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
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Edge(5),
            legs: &legs,
        };
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
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(0),
            legs: &legs,
        };
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
        let incompatible = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(0),
            legs: &legs,
        };
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
            let instantiation = NumeratorInstantiation {
                owner: NumeratorOwner::Vertex(0),
                legs: &legs,
            };
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
        let closed = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(3),
            legs: &[],
        };
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
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(2),
            legs: &legs,
        };
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
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(9),
            legs: &legs,
        };
        let numerator = generator
            .interaction_numerator(model.vertex_rule_id("V_78").unwrap(), &instantiation)
            .unwrap();

        let numerator = NumeratorInterface::export(&numerator).0.to_plain_string();
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
        let signatures = InteractionSignatures::new(&model, &options).unwrap();
        let rules: BTreeSet<_> = signatures.0.values().flatten().copied().collect();
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
            combined.push_str(
                &NumeratorInterface::export(diagram.numerator())
                    .0
                    .to_plain_string(),
            );
            for (_, vertex) in diagram.vertices() {
                combined.push_str(
                    &NumeratorInterface::export(&vertex.numerator)
                        .0
                        .to_plain_string(),
                );
            }
            for (_, _, edge) in diagram.edges() {
                combined.push_str(
                    &NumeratorInterface::export(&edge.numerator)
                        .0
                        .to_plain_string(),
                );
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
        assert_eq!(diagram.numerator().expression(), &Atom::one());
        let routed = diagram.loop_momentum_basis().clone();
        let id = diagram.id();
        let diagram = diagram.apply_feynman_rules().unwrap();
        assert_eq!(diagram.numerator().expression(), &test_atom("GC1^2"));
        assert_eq!(diagram.loop_momentum_basis(), &routed);
        assert_eq!(diagram.id(), id);
        diagram.validate().unwrap();
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
            let topology = diagram
                .clone()
                .with_numerator(SymbolicTensor::infer(Atom::one()).unwrap())
                .unwrap();
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
                EdgeColor {
                    particle: particle(&generator.model, 1),
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
                    particle: particle(&generator.model, 22),
                    direction: None,
                },
            )
            .unwrap();

        let loop_count = generator.closed_fermion_loop_count(&graph, false).unwrap();
        assert_eq!(loop_count, 1);
        let diagram = generator
            .to_diagram(
                "loop".to_owned(),
                GenerationType::Amplitude,
                ColoredTopology {
                    graph,
                    symmetry: 3,
                    multiplicity: 2,
                    cut_partitions: Vec::new(),
                    topology_threshold_candidates: Vec::new(),
                },
                FermionSigns {
                    internal_loop_negative: loop_count % 2 == 1,
                    ..FermionSigns::default()
                },
                &GenerationOptions::default(),
            )
            .unwrap();
        assert_eq!(
            diagram.overall_factor(),
            &factor_atom("InternalFermionLoopSign(-1)*CouplingsMultiplicity(2)/AutG(3)")
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
                    .add_edge(
                        left,
                        right,
                        true,
                        EdgeColor {
                            particle: particle(&model, pdg),
                            direction: None,
                        },
                    )
                    .unwrap();
            }
            signs.push(
                process
                    .normalize_fermion_flows(&graph, &model)
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
        // This artificial spin/statistics model no longer describes the SM's
        // physical vector/Goldstone/ghost quartets used for covariant cuts.
        definition["covariant_cut_multiplets"] = serde_json::json!({});
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
                    generator.closed_fermion_loop_count(&graph, true).unwrap(),
                    expected,
                    "{name}"
                );
                assert_eq!(
                    generator.closed_fermion_loop_count(&graph, false).unwrap(),
                    expected_fermions,
                    "fermion filter: {name}"
                );
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
        graph
            .add_edge(incoming, outgoing, true, edge(&generator.model, 9000005))
            .unwrap();
        assert_eq!(
            generator.closed_fermion_loop_count(&graph, true).unwrap(),
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
                generator.closed_fermion_loop_count(&graph, true).unwrap(),
                loop_count,
            );
            assert_eq!(
                generator.closed_fermion_loop_count(&graph, false).unwrap(),
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
                "positive".to_owned(),
                GenerationType::Amplitude,
                ColoredTopology {
                    graph: Graph::new(),
                    symmetry: 1,
                    multiplicity: 1,
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
                .add_edge(
                    center,
                    endpoint,
                    true,
                    EdgeColor {
                        particle: particle(&generator.model, 1),
                        direction: None,
                    },
                )
                .unwrap();
        }

        assert!(matches!(
            generator.closed_fermion_loop_count(&graph, false),
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
        let canonical_only = GenerationOptions::default().numerator_grouping(
            NumeratorGrouping::Identical(GraphGroupingOptions {
                number_of_numerical_samples: 0,
                check_canonical_numerator: true,
                ..GraphGroupingOptions::default()
            }),
        );
        generator
            .validate_options(GenerationType::Amplitude, &canonical_only)
            .unwrap();
    }

    #[test]
    fn public_limits_report_arithmetic_overflow() {
        let model = fermion_model();
        let generator = Generator::new(model.clone());
        assert!(matches!(
            generator.generate(
                &photon_process(),
                &GenerationOptions::default().max_vertices(usize::MAX),
                GenerationType::Amplitude
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
                    .cut_partitions(&graph, &model, &options)
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
    fn external_canonical_classes_keep_global_left_right_symmetry_out_of_raw_generation() {
        let model = fermion_model();
        let classes = |process: Process, generation_type, options: GenerationOptions| {
            ResolvedProcess::new(&model, &process, generation_type, &options)
                .unwrap()
                .external_edges(&model)
                .unwrap()
                .into_iter()
                .map(|(node, _)| node.external_class.unwrap())
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
                ExternalCanonicalClass::Exact {
                    state: ExternalState::Incoming,
                    index: 0,
                },
                ExternalCanonicalClass::Exact {
                    state: ExternalState::Incoming,
                    index: 1,
                },
                ExternalCanonicalClass::WithinState(ExternalState::Incoming),
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
                ExternalCanonicalClass::WithinState(ExternalState::Incoming),
                ExternalCanonicalClass::WithinState(ExternalState::Incoming),
                ExternalCanonicalClass::WithinState(ExternalState::Incoming),
            ]
        );

        let cross_section = Process::new([1_i64, -1], [22_i64]);
        assert_eq!(
            classes(
                cross_section,
                GenerationType::CrossSection,
                GenerationOptions::default().symmetrize_initial(true)
            ),
            vec![
                ExternalCanonicalClass::Exact {
                    state: ExternalState::Incoming,
                    index: 0
                },
                ExternalCanonicalClass::Exact {
                    state: ExternalState::Incoming,
                    index: 1
                },
                ExternalCanonicalClass::Exact {
                    state: ExternalState::Outgoing,
                    index: 2,
                },
                ExternalCanonicalClass::Exact {
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
                ExternalCanonicalClass::Exact {
                    state: ExternalState::Incoming,
                    index: 0,
                },
                ExternalCanonicalClass::Exact {
                    state: ExternalState::Incoming,
                    index: 1,
                },
                ExternalCanonicalClass::Exact {
                    state: ExternalState::Outgoing,
                    index: 2,
                },
                ExternalCanonicalClass::Exact {
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
                EdgeColor {
                    particle: particle(&generator.model, 1),
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
                    particle: particle(&generator.model, 22),
                    direction: None,
                },
            )
            .unwrap();

        let diagram = generator
            .to_diagram(
                "localized".to_owned(),
                GenerationType::Amplitude,
                ColoredTopology {
                    graph,
                    symmetry: 1,
                    multiplicity: 1,
                    cut_partitions: Vec::new(),
                    topology_threshold_candidates: Vec::new(),
                },
                FermionSigns {
                    internal_loop_negative: true,
                    ..FermionSigns::default()
                },
                &GenerationOptions::default(),
            )
            .unwrap();
        let internal_numerator = &diagram
            .edges()
            .find(|(_, _, edge)| edge.external.is_none())
            .unwrap()
            .2
            .numerator;
        let internal_numerator = NumeratorInterface::export(internal_numerator)
            .0
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
        assert_eq!(external_numerator.expression(), &Atom::one());

        let vertex_numerator = &diagram
            .vertex(feynkit_graph::VertexId(interaction))
            .unwrap()
            .numerator;
        let vertex_numerator = NumeratorInterface::export(vertex_numerator)
            .0
            .to_plain_string();
        assert!(vertex_numerator.contains("gammalooprs::hedge"));
        assert!(vertex_numerator.contains("gammalooprs::hedge"));
    }
}
