use feynkit_graph::symbols;
use idenso::tensor::composition;
use idenso::tensor::inference::InterfaceInference;
use idenso::{CookMode, CookSettings};
use spenso::structure::abstract_index::AIND_SYMBOLS;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialIndex, PartialStructureExt};
use spenso::structure::representation::LibraryRep;
use spenso::structure::slot::{IsAbstractSlot, ParseableAind, SlotMatch, SlotMatcher};
use symbolica::atom::UserData;
// The case table of the generation oracle tests, shared by the unit tests that run its workloads.
#[cfg(test)]
#[path = "../tests/generation_oracle/cases.rs"]
#[allow(dead_code)]
mod cases;
mod enumeration;
mod graph;
mod interaction_stars;
#[cfg(test)]
mod projector_equivalence;
#[cfg(test)]
mod search_census;
#[cfg(test)]
mod signed_routing;
mod special_topologies;
mod staged;
pub use staged::{BuiltDiagramSet, Counts, DiagramSet, LoopCounts, PreparedDiagram};

#[cfg(test)]
use std::sync::atomic::{AtomicUsize, Ordering};

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt;
use std::hash::Hash;
use std::ops::ControlFlow;
use std::sync::{Arc, LazyLock, Mutex, RwLock};

use feynkit_graph::{
    DiagramCut, DiagramCutSide, DiagramEdge, DiagramEndpoint, DiagramError, DiagramHalfEdge,
    DiagramId, DiagramThresholdCandidate, DiagramVertex, DotDiagram, EdgeId, ExternalLeg,
    ExternalState, FeynmanDiagram, FeynmanDiagramBuilder, NumeratorInterface, VertexId, VertexSlot,
};
#[cfg(test)]
use feynkit_model::VertexRule;
use feynkit_model::{Model, ModelError, Particle, ParticleId, VertexRuleId};
use idenso::color::CS;
use idenso::dirac::AGS;
use idenso::epsilon::EPSILON_SYMBOL;
use idenso::representations::Bispinor;
use idenso::tensor::SymbolicTensor;
use linnet::half_edge::canon::{CanonicalForm, KeyedEdge};
use linnet::half_edge::generate::{
    Diagram as CompactDiagram, Plan, Run as EngineRun, SkeletonSnapshot,
};
use linnet::half_edge::involution::{Flow, Orientation};
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike};
use linnet::half_edge::{HedgeGraph, NoData, NodeIndex};
use rayon::ThreadPoolBuildError;
#[cfg(not(target_arch = "wasm32"))]
use rayon::ThreadPoolBuilder;
use serde::{Deserialize, Serialize};
use spenso::network::library::symbolic::ETS;
use spenso::structure::partial::PartialStructure;
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use symbolica::symbol;
use thiserror::Error;

use feynkit_amplitude::ColorRepresentation;

use crate::options::ResolvedRestrictions;
use crate::{
    GenerationType, Grouping, GroupingError, Naming, NumeratorOptions, ParticleSelector, Process,
    ProcessError, Request, Run, Selection, SelectorError, Stage, momentum_symbol,
};
pub use enumeration::{
    ExternalLabel, GenerationDiagramView, GenerationDressingView, GenerationSkeletonView,
};
use enumeration::{
    ParticleStub, ProcessExternals, Pruner, VertexSignatures, closed_loops, incoming_stub,
};
pub use graph::{Endpoint as PreparedEndpoint, ExternalNode, GeneratedEdge as PreparedEdge};
use graph::{Endpoint, GeneratedEdge, GeneratedGraph, PlainNodes, PlainVertex};
use interaction_stars::{InteractionStar, InteractionStars};
use num_bigint::{BigInt, BigUint};
use num_rational::BigRational;
#[cfg(test)]
use signed_routing::{DiagramWeight, SignedRoutedDiagram};

#[cfg(test)]
#[path = "main_generation_regressions.rs"]
mod main_generation_regressions;

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
    #[error("graph automorphism factor cannot be represented by u128")]
    SymmetryFactor,
    #[error(transparent)]
    Enumeration(#[from] linnet::half_edge::generate::GenerationError),
    #[error("invalid generation request: {0}")]
    InvalidRequest(String),
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
    /// Index into [`BuiltDiagramSet::diagrams`] after zero-numerator removal.
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

impl DiagramGroup {
    /// The group of one diagram, its own master.
    pub(crate) fn singleton(
        diagram: usize,
        source_diagram: usize,
        source_identity: Option<SourceIdentity>,
        overall_factor: Atom,
    ) -> Self {
        Self {
            master: diagram,
            members: vec![GroupMember {
                source_diagram,
                source_identity,
                diagram,
                ratio: Atom::one(),
                overall_factor,
            }],
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenerationReport {
    /// Complete untyped topologies accepted during enumeration. Loop orders that the coupling
    /// orders rule out (no multiset of vertex rules has the degrees and orders they need), and
    /// requests with external charges or parities that no vertex rule balances or with an empty
    /// coupling-order or fermion-loop range, are not enumerated and add none.
    pub topology_count: BigUint,
    /// Enumerated diagrams: topologies with fields and vertex rules, before the cut and
    /// charge-conjugation filters that need their canonical frame.
    pub interaction_assignment_count: BigUint,
    /// Number of diagrams retained after selection and any requested numerator filtering.
    pub retained_count: BigUint,
    pub zero_numerator_count: BigUint,
    pub completed: bool,
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
    options: &Request,
) -> Result<Option<UnresolvedCutContent>, GenerationError> {
    let required = &options.orders.perturbative;
    if required.is_empty() {
        return Ok(None);
    }
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

impl<'legs> NumeratorInstantiation<'legs> {
    fn new(owner: NumeratorOwner, legs: &'legs [NumeratorHalfEdge]) -> Self {
        Self { owner, legs }
    }

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

/// Resolved process inputs and their one shared search builder. The search census and
/// projector oracle use the same pruning and budget construction as ordinary generation.
struct PreparedSearch {
    request: Request,
    restrictions: ResolvedRestrictions,
    resolved: ResolvedProcess,
    signatures: VertexSignatures,
    externals: ProcessExternals,
    plan: Option<Plan<ParticleId, ExternalLabel>>,
    report: linnet::half_edge::generate::PlanReport,
}

#[derive(Clone)]
pub(crate) struct Generator {
    model: Arc<Model>,
    cache: Arc<GenerationCache>,
    numerators: Arc<NumeratorCache>,
}

/// Reuse model-local interaction stars, canonical routing and exact symbolic factors
/// for this generator's lifetime. Per-diagram identities are bound at each use.
struct GenerationCache {
    stars: InteractionStars,
    routing: Arc<feynkit_graph::routing::CanonicalRoutingCache>,
    factors: RwLock<HashMap<(u128, u64, FermionSigns), Atom>>,
    /// Automatic projectors by the finalized position of every external state. Finalization
    /// numbers the external legs first, so all diagrams of an amplitude share one entry. A cross
    /// section adds the sewn final-state copies after the internal edges, so its entries vary
    /// with the internal edge count.
    projectors: RwLock<HashMap<Vec<ExternalStatePosition>, SymbolicTensor<PartialStructure>>>,
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

/// Localized factors depend on the immutable model and their index assignment,
/// not on the rest of the diagram. Keep them for this generator's lifetime.
#[derive(Default)]
struct NumeratorCache {
    interactions: RwLock<HashMap<InteractionNumeratorKey, SymbolicTensor<PartialStructure>>>,
    propagators: RwLock<HashMap<PropagatorNumeratorKey, SymbolicTensor<PartialStructure>>>,
    interaction_templates: RwLock<HashMap<InteractionNumeratorKey, Arc<NumeratorTemplate>>>,
    propagator_templates: RwLock<HashMap<ParticleId, Arc<NumeratorTemplate>>>,
}

#[derive(Clone, Copy)]
enum InteractionMomenta {
    Unreduced,
    Conserved,
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
        let Self::Admitted { tensor, .. } = self else {
            unreachable!("dynamic rules are constructed at their actual owner")
        };
        NumeratorInterface::translate(tensor, |atom| instantiation.bind_template(atom))
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

    /// Construct the external wavefunctions in the finalized half-edge frame: the
    /// projector that generation and [`Self::apply_feynman_rules`] build for an omitted
    /// one. Callers with an explicit projector retain that supplied expression.
    fn external_state_projector(&self)
    -> Result<SymbolicTensor<PartialStructure>, GenerationError>;
}

impl FeynmanDiagramRulesExt for FeynmanDiagram {
    fn apply_feynman_rules(&self) -> Result<FeynmanDiagram, GenerationError> {
        // `map_data` returns unvalidated diagrams: check their slots before they
        // index the vertex rules. The check is cached for validated diagrams.
        self.validate()?;
        Generator::new(self.model_arc()).construct_numerators(self.clone(), |_| true)
    }

    fn external_state_projector(
        &self,
    ) -> Result<SymbolicTensor<PartialStructure>, GenerationError> {
        self.validate()?;
        Generator::new(self.model_arc()).finalized_external_state_projector(self)
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

    fn external_state_projector(
        &self,
    ) -> Result<SymbolicTensor<PartialStructure>, GenerationError> {
        self.diagram.external_state_projector()
    }
}

impl Generator {
    pub fn new(model: impl Into<Arc<Model>>) -> Self {
        let model = model.into();
        let cache = Arc::new(GenerationCache {
            stars: InteractionStars::new(&model),
            routing: Arc::default(),
            factors: RwLock::default(),
            projectors: RwLock::default(),
        });
        Self {
            model,
            cache,
            numerators: Arc::default(),
        }
    }

    fn search(
        &self,
        request: &Request,
        observe_skeleton: impl Fn(&GenerationSkeletonView<'_>, bool) + Send + Sync + 'static,
        observe_dressing: impl Fn(&GenerationDressingView<'_>, bool) + Send + Sync + 'static,
    ) -> Result<PreparedSearch, GenerationError> {
        request.validate_in(&self.model)?;
        let restrictions = request.restrictions.resolve(&self.model)?;
        let expanded = request.process.clone().with_final_state_alternatives(
            request
                .process
                .covariant_cut_states(&self.model, request.kind)?,
        )?;
        let resolved = ResolvedProcess::new(&self.model, &expanded, request.kind, request)?;
        let signatures = VertexSignatures::new(&self.model, request)?;
        let externals = resolved.externals(&self.model)?;
        let pruner = Arc::new(Pruner::new(&self.model, request));
        let skeleton_pruner = Arc::clone(&pruner);
        let observe = request.deferred_hooks.is_none();
        let skeleton_hook = request.prune_skeletons.clone();
        let dressing_hook = request.prune_dressings.clone();
        let planned = signatures.plan(
            &self.model,
            &externals,
            request,
            move |view| {
                let keep = skeleton_pruner.keeps_skeleton(view)
                    && skeleton_hook.as_ref().is_none_or(|f| f(view));
                if observe {
                    observe_skeleton(view, keep)
                }
                keep
            },
            move |view| {
                let keep =
                    pruner.keeps_dressing(view) && dressing_hook.as_ref().is_none_or(|f| f(view));
                if observe {
                    observe_dressing(view, keep)
                }
                keep
            },
        );
        let (plan, report) = match planned {
            Ok(plan) => {
                let report = plan.explain().clone();
                (Some(plan), report)
            }
            Err(GenerationError::Enumeration(
                linnet::half_edge::generate::GenerationError::Infeasible(report),
            )) => (None, report),
            Err(error) => return Err(error),
        };
        Ok(PreparedSearch {
            request: request.clone(),
            restrictions,
            resolved,
            signatures,
            externals,
            plan,
            report,
        })
    }

    /// Canonicalize a generated diagram, apply the cut and interaction filters that need its
    /// canonical frame, and normalize its fermion flows.
    fn diagram_input(
        &self,
        resolved: &ResolvedProcess,
        graph: GeneratedGraph,
        symmetry: u128,
        multiplicity: u64,
        options: &Request,
        grouping: &Grouping,
    ) -> Result<Option<DiagramInput>, GenerationError> {
        let mut topology = ColoredTopology {
            graph: resolved.canonicalize_numerator_graph(&graph, &self.model, grouping)?,
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
        topology.graph = self.repair_cp_interaction_rules(topology.graph)?;
        if !self.passes_vertex_selection(&topology.graph, options)? {
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
        let (comparison, reversed_edges) = if matches!(*grouping, Grouping::None) {
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

    /// Finalize the comparison and representative diagrams of a prepared input. Equal frames
    /// with equal requested prefactors and projectors share one finalized owner.
    fn diagram_pair(
        &self,
        name: Option<String>,
        generation_type: GenerationType,
        input: DiagramInput,
        options: (&NumeratorOptions, &NumeratorOptions),
    ) -> Result<DiagramPair, GenerationError> {
        let DiagramInput {
            comparison,
            representative,
            signs,
            reversed_edges,
            multiplicity,
        } = input;
        let identical = comparison == representative
            && options.0.prefactor == options.1.prefactor
            && options.0.projector == options.1.projector;
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
            options.0,
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
            options.1,
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

    /// The finalized physical geometry of a prepared input with its exact signed weight,
    /// without local Feynman rules or an external projector.
    #[cfg(test)]
    fn routed_diagram(
        &self,
        name: Option<String>,
        generation_type: GenerationType,
        input: DiagramInput,
    ) -> Result<SignedRoutedDiagram, GenerationError> {
        let weight = DiagramWeight::new(input.signs, input.multiplicity);
        let builder = FeynmanDiagramBuilder::new(Arc::clone(&self.model), name);
        let (geometry, _) =
            self.assemble_diagram(generation_type, input.representative, builder)?;
        Ok(SignedRoutedDiagram::new(geometry, weight))
    }

    /// Recheck allowed rules after charge-conjugation repair.
    fn passes_vertex_selection(
        &self,
        graph: &GeneratedGraph,
        request: &Request,
    ) -> Result<bool, GenerationError> {
        let restrictions = request.restrictions.resolve(&self.model)?;
        Ok(graph.vertices().all(|rule| restrictions.vertex(rule)))
    }

    fn retain_cut_amplitude_partitions(
        &self,
        topology: &mut ColoredTopology,
        request: &Request,
    ) -> Result<bool, GenerationError> {
        let required = &request.cuts.amplitude;
        let mut retained = Vec::with_capacity(topology.cut_partitions.len());
        for partition in topology.cut_partitions.drain(..) {
            let total_loops = partition
                .left
                .loop_count
                .checked_add(partition.right.loop_count)
                .ok_or(GenerationError::ArithmeticOverflow("cut-amplitude loops"))?;
            let passes = required.loops.contains(total_loops)
                && [&partition.left, &partition.right].into_iter().all(|side| {
                    required.orders.plain.iter().all(|(name, range)| {
                        range.contains(side.coupling_orders.get(name).copied().unwrap_or(0))
                    })
                });
            if passes {
                retained.push(partition)
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
        options: &NumeratorOptions,
    ) -> Result<FeynmanDiagram, GenerationError> {
        let overall_factor = memoised(
            &self.cache.factors,
            (topology.symmetry, multiplicity, fermion_signs),
            |&(symmetry, multiplicity, signs)| Ok(signs.factor(symmetry, multiplicity)),
        )?;
        let mut builder = FeynmanDiagramBuilder::new(Arc::clone(&self.model), name)
            .overall_factor(overall_factor);
        if !options.prefactor.is_one() {
            builder = builder
                .numerator_prefactor(NumeratorInterface::admit(options.prefactor.clone(), None)?);
        }
        let (diagram, external_states) =
            self.assemble_diagram(generation_type, topology, builder)?;
        let positions = external_states
            .into_iter()
            .map(|(external, connection)| {
                ExternalStatePosition::finalized(external, connection, &diagram)
            })
            .collect::<Vec<_>>();
        let projector = match &options.projector {
            Some(supplied) => Self::supplied_external_state_projector(supplied, &positions)?,
            None => memoised(&self.cache.projectors, positions, |positions| {
                self.external_state_projector(positions)
            })?,
        };
        self.construct_numerators(diagram.with_projector(projector), |_| true)
    }

    /// Materialize and route the same physical geometry before any Feynman rules.
    fn assemble_diagram(
        &self,
        generation_type: GenerationType,
        topology: ColoredTopology,
        builder: FeynmanDiagramBuilder,
    ) -> Result<(FeynmanDiagram, Vec<(ExternalNode, usize)>), GenerationError> {
        let graph = &topology.graph;
        let incoming_count = graph
            .externals()
            .filter(|(_, external, _)| external.state == ExternalState::Incoming)
            .count();
        let mut builder = builder
            .canonical_loop_momentum_basis()
            .canonical_routing_cache(Arc::clone(&self.cache.routing))
            .symmetry_factor(topology.symmetry)
            .cuts(topology.cut_partitions.clone())
            .topology_threshold_candidates(topology.topology_threshold_candidates.clone());
        // Every external state with the connection of the leg that carries it.
        let mut external_states = Vec::new();
        let mut endpoint_slots = vec![[VertexSlot(0); 2]; graph.n_edges()];
        let mut star = InteractionStar::default();
        for (vertex, rule_id) in graph.vertices().enumerate() {
            let rule = self.model.vertex_rule_by_id(rule_id)?;
            self.cache
                .stars
                .read(graph, vertex, &self.model, &mut star)?;
            self.cache.stars.rule(rule_id).bind(
                &star,
                vertex,
                rule,
                &self.model,
                |slot, leg| {
                    endpoint_slots[leg.edge][usize::from(leg.flow.is_sink())] = VertexSlot(slot);
                },
            )?;
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
        Ok((builder.build()?, external_states))
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
            let ports = self.cache.stars.rule(rule_id);
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
                let (spin, color) = ports.spin_color(slot.0);
                legs[slot.0] = Some(NumeratorHalfEdge {
                    edge: edge.0,
                    hedge: hedge.0,
                    flow,
                    spin,
                    color,
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
    ) -> Result<SymbolicTensor<PartialStructure>, GenerationError> {
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
        Ok(NumeratorInterface::admit(projector, None)?)
    }

    /// A requested projector in finalized indices. Process-level expressions name external
    /// legs by their process index: an edge index or carrier names the edge of that leg and a
    /// half-edge index the half-edge where it attaches. Native local UFO slots are unchanged.
    fn supplied_external_state_projector(
        supplied: &Atom,
        positions: &[ExternalStatePosition],
    ) -> Result<SymbolicTensor<PartialStructure>, GenerationError> {
        let positions = positions
            .iter()
            .map(|position| (position.external.index, (position.edge, position.hedge)))
            .collect::<BTreeMap<_, _>>();
        let hedge_head = feynkit_graph::symbols::hedge_index();
        let edge_head = feynkit_graph::symbols::edge_index();
        let carriers = [
            momentum_symbol(),
            feynkit_graph::symbols::u(),
            feynkit_graph::symbols::ubar(),
            feynkit_graph::symbols::v(),
            feynkit_graph::symbols::vbar(),
            feynkit_graph::symbols::epsilon(),
            feynkit_graph::symbols::epsilonbar(),
        ];
        let mut projector = supplied.clone();
        // Resolve index owners before their carriers: a replaced carrier is not visited again.
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
                let Some(&(edge, hedge)) = positions.get(&index) else {
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

    /// The external-state projector of a finalized diagram. An incoming state is
    /// indexed at the target of its edge and an outgoing one at the source; a sewn
    /// cross-section initial state carries both copies.
    fn finalized_external_state_projector(
        &self,
        diagram: &FeynmanDiagram,
    ) -> Result<SymbolicTensor<PartialStructure>, GenerationError> {
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
        Ok(NumeratorInterface::admit(projector, None)?)
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

    #[cfg(test)]
    fn interaction_half_edges(
        &self,
        graph: &GeneratedGraph,
        vertex: usize,
        rule: &VertexRule,
    ) -> Result<Vec<NumeratorHalfEdge>, GenerationError> {
        let rule_id = self.model.vertex_rule_id(&rule.name)?;
        let mut star = InteractionStar::default();
        self.cache
            .stars
            .read(graph, vertex, &self.model, &mut star)?;
        let mut legs = Vec::with_capacity(rule.particles.len());
        self.cache
            .stars
            .rule(rule_id)
            .bind(&star, vertex, rule, &self.model, |_, leg| legs.push(leg))?;
        Ok(legs)
    }

    fn repair_cp_interaction_rules(
        &self,
        graph: GeneratedGraph,
    ) -> Result<GeneratedGraph, GenerationError> {
        let mut replacements = BTreeMap::new();
        let mut star = InteractionStar::default();
        for (vertex, rule_id) in graph.vertices().enumerate() {
            self.model.vertex_rule_by_id(rule_id)?;
            self.cache
                .stars
                .read(&graph, vertex, &self.model, &mut star)?;
            let repaired = self
                .cache
                .stars
                .repair_rule(rule_id, vertex, &star, &self.model)?;
            if repaired != rule_id {
                replacements.insert(vertex, repaired);
            }
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
            // Keep one actual-owner admission for callback and compact rules.
            // Reconstructing them would replay normalization or allocate new
            // dummy identities. Conservation shortening belongs to the retained
            // intrinsic template, not a second lowering of this dynamic rule.
            NumeratorTemplate::Dynamic => NumeratorInterface::admit(
                self.interaction_expression(rule_id, instantiation, InteractionMomenta::Unreduced)?,
                None,
            )?,
            NumeratorTemplate::Admitted { .. } => template.instantiate(instantiation)?,
        }
        // Raw rules declare no axis order. Named interfaces use the same
        // canonical convention after actual-owner admission and template
        // transport; leaf argument order continues to specify the tensor.
        .canonicalize_named_interface();
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
                let expression = self.interaction_expression(
                    rule_id,
                    &instantiation,
                    InteractionMomenta::Unreduced,
                )?;
                // Admit every raw summand and retain the full label universe before
                // conservation can erase the payload. Classification also belongs
                // to this unreduced boundary, just as it does before edge aliasing.
                let template =
                    NumeratorTemplate::admitted(NumeratorInterface::admit(expression, None)?)?;
                match template {
                    NumeratorTemplate::Admitted { tensor, labels } => NumeratorTemplate::Admitted {
                        tensor: Box::new(self.reduce_interaction_momenta(
                            &tensor,
                            rule_id,
                            &instantiation,
                        )?),
                        labels,
                    },
                    NumeratorTemplate::Dynamic => NumeratorTemplate::Dynamic,
                }
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

    /// Rewrite only after the raw interaction has established its interface and
    /// construction history. Cooking changes notation, not tensor admission;
    /// the existing checked rewrite retains hidden indices and typed zeros.
    fn reduce_interaction_momenta(
        &self,
        tensor: &SymbolicTensor<PartialStructure>,
        rule_id: VertexRuleId,
        instantiation: &NumeratorInstantiation<'_>,
    ) -> Result<SymbolicTensor<PartialStructure>, GenerationError> {
        let expression =
            self.interaction_expression(rule_id, instantiation, InteractionMomenta::Conserved)?;
        let expression = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true)
            .try_cook(expression.as_view())
            .map_err(|source| DiagramError::TensorCooking {
                operation: "encoding the momentum-conserved interaction",
                source: Box::new(source),
            })?;
        Ok(tensor
            .with_checked_expression(expression)
            .map_err(DiagramError::from)?)
    }

    fn interaction_expression(
        &self,
        rule_id: VertexRuleId,
        instantiation: &NumeratorInstantiation<'_>,
        momenta: InteractionMomenta,
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
                    let lorentz = self.model.lorentz_structure_by_id(*lorentz_name)?;
                    let conserved;
                    let lorentz = match momenta {
                        InteractionMomenta::Unreduced => &lorentz.structure,
                        InteractionMomenta::Conserved => {
                            conserved = lorentz.reduced_by_momentum_conservation();
                            &conserved
                        }
                    };
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
            NumeratorTemplate::Admitted { .. } => template.instantiate(&instantiation)?,
        }
        .canonicalize_named_interface();
        self.numerators
            .propagators
            .write()
            .expect("numerator cache was poisoned; its retained tensors cannot be reused")
            .insert(key, numerator.clone());
        Ok(numerator)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct ColoredTopology {
    graph: GeneratedGraph,
    symmetry: u128,
    cut_partitions: Vec<DiagramCut>,
    topology_threshold_candidates: Vec<DiagramThresholdCandidate>,
}

/// A canonicalized diagram ready for finalization: the frame used for numerator comparisons,
/// the physical representative with normalized fermion flows, and the edges reversed between
/// them. Both frames reframe the same generated assignment, so [`GeneratedGraph::frame_map`]
/// relates either to any other preparation of it.
#[derive(Clone)]
struct DiagramInput {
    comparison: ColoredTopology,
    representative: ColoredTopology,
    signs: FermionSigns,
    reversed_edges: BTreeSet<EdgeId>,
    /// The number of exactly labelled generated diagrams the representative stands for.
    multiplicity: u64,
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

    /// Retain `cuts`, given in the physical representative's frame, in both frames.
    fn with_cuts(self, cuts: Vec<DiagramCut>) -> Result<Self, DiagramError> {
        match self {
            Self::Shared(diagram) => Ok(Self::Shared(diagram.with_cuts(cuts)?)),
            Self::Distinct {
                comparison,
                representative,
                reversed_edges,
            } => {
                let representative = (*representative).with_cuts(cuts)?;
                let partitions = transport_cut_metadata(representative.cuts(), &reversed_edges)
                    .into_iter()
                    .map(|cut| (cut.left.half_edges, cut.right.half_edges))
                    .collect();
                Ok(Self::Distinct {
                    comparison: comparison.with_cut_partitions(partitions)?,
                    representative: Box::new(representative),
                    reversed_edges,
                })
            }
        }
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

/// The bookkeeping heads of the overall factors, created together in this order on first use.
/// Symbolica orders the factors of a product by symbol creation, so a head first created by
/// whichever generation worker needs it first would make printed factors depend on the
/// thread schedule.
pub(crate) struct FactorHeads {
    aut_g: Symbol,
    numerator_independent_grouping: Symbol,
    internal_fermion_loop_sign: Symbol,
    external_fermion_ordering_sign: Symbol,
    antifermion_spin_sum_sign: Symbol,
    pub(crate) numerator_dependent_grouping: Symbol,
}

pub(crate) static FACTOR_HEADS: LazyLock<FactorHeads> = LazyLock::new(|| FactorHeads {
    aut_g: symbol!("feynkit_generator_factor::AutG"),
    numerator_independent_grouping: symbol!(
        "feynkit_generator_factor::NumeratorIndependentSymmetryGrouping"
    ),
    internal_fermion_loop_sign: symbol!("feynkit_generator_factor::InternalFermionLoopSign"),
    external_fermion_ordering_sign: symbol!(
        "feynkit_generator_factor::ExternalFermionOrderingSign"
    ),
    antifermion_spin_sum_sign: symbol!("feynkit_generator_factor::AntiFermionSpinSumSign"),
    numerator_dependent_grouping: symbol!("feynkit_generator_factor::NumeratorDependentGrouping"),
});

impl FermionSigns {
    /// The diagram multiplier depends only on its exact weight and physical
    /// signs, so decorated graphs with the same values share its construction.
    fn factor(self, symmetry: u128, multiplicity: u64) -> Atom {
        let heads = &*FACTOR_HEADS;
        // The automorphism count of the rule-coloured graph is a denominator
        // in the overall multiplier.
        let mut factors = vec![heads.aut_g.call(symmetry).pow(-1)];
        if multiplicity > 1 {
            factors.push(heads.numerator_independent_grouping.call(multiplicity));
        }
        if self.internal_loop_negative {
            factors.push(heads.internal_fermion_loop_sign.call(-1));
        }
        if self.include_external_ordering {
            let sign = if self.external_ordering_negative {
                -1
            } else {
                1
            };
            factors.push(heads.external_fermion_ordering_sign.call(sign));
        }
        if self.include_antifermion_spin_sum {
            let sign = if self.antifermion_spin_sum_negative {
                -1
            } else {
                1
            };
            factors.push(heads.antifermion_spin_sum_sign.call(sign));
        }
        Atom::mul_many(factors)
    }
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
        options: &Request,
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
            symmetrize_initial: options.symmetrization.initial,
            symmetrize_final: options.symmetrization.final_state,
            symmetrize_left_right: options.cuts.symmetrize_left_right,
            symmetrize_external_fermions: options.symmetrization.external_fermions,
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
        grouping: &Grouping,
    ) -> Result<GeneratedGraph, GenerationError> {
        let grouping_options = match grouping {
            Grouping::Identical(options)
            | Grouping::UpToSign(options)
            | Grouping::UpToScalar(options) => Some(options),
            Grouping::None => None,
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
        if include_left_right {
            self.canonicalize_numerator_graph_pass(
                &selected,
                model,
                exact_internal_particles,
                false,
                true,
            )
        } else {
            Ok(selected)
        }
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
        // Internal fermion direction is omitted from the reduced skeleton but restored below.
        // Ghost vertices differentiate one propagator end, so a ghost edge keeps its
        // ghost-number arrow: reversing storage while conjugating the ghost describes the
        // same directed edge. The two global left/right choices share incidence and edge
        // colors. Prepare that skeleton once and vary only its external half-edge keys when
        // canonizing it.
        let edges = graph
            .edges()
            .map(|edge| {
                let color = self
                    .numerator_colors
                    .edge(&edge, model, exact_internal_particles)?;
                let particle = model.particle_by_id(edge.particle)?;
                let labels = match (edge.source, edge.sink) {
                    (Endpoint::External(external), _) | (_, Endpoint::External(external)) => {
                        let normal = external_color(external)?;
                        let swapped = if include_left_right_remapping {
                            external_color(partner(external))?
                        } else {
                            normal
                        };
                        Some([normal, swapped])
                    }
                    _ => None,
                };
                Ok((edge, color, labels, particle))
            })
            .collect::<Result<Vec<_>, GenerationError>>()?;
        // The reduced skeleton of one global left/right choice, in edge order.
        let reduced_edges = |swapped: bool| {
            edges.iter().map(
                move |&(edge, key, labels, particle)| match (edge.source, edge.sink) {
                    (Endpoint::Vertex(source), Endpoint::Vertex(sink)) => {
                        let ghost = particle.is_ghost().then(|| particle.is_antiparticle());
                        let ((source, sink), orientation) = match ghost {
                            Some(true) => ((sink, source), Orientation::Default),
                            Some(false) => ((source, sink), Orientation::Default),
                            None => ((source, sink), Orientation::Undirected),
                        };
                        KeyedEdge::Paired {
                            source: (NodeIndex(source), None),
                            sink: (NodeIndex(sink), None),
                            key,
                            orientation,
                        }
                    }
                    (Endpoint::External(_), Endpoint::Vertex(vertex))
                    | (Endpoint::Vertex(vertex), Endpoint::External(_)) => KeyedEdge::Dangling {
                        node: NodeIndex(vertex),
                        hedge_key: labels.map(|labels| labels[usize::from(swapped)]),
                        key,
                        orientation: Orientation::Default,
                        flow: Flow::Sink,
                    },
                    (Endpoint::External(_), Endpoint::External(_)) => {
                        unreachable!("generated edges touch an interaction vertex")
                    }
                },
            )
        };
        // The physical graph of one global left/right choice, in edge order: fields in their
        // physical direction, and legs keyed by the external index they show.
        let physical_edges = |swapped: bool| {
            edges.iter().map(move |&(edge, _, _, particle)| {
                let reverse = edge.directed && particle.is_antiparticle();
                let key = if reverse {
                    particle.antiparticle
                } else {
                    edge.particle
                };
                let orientation = if !edge.directed {
                    Orientation::Undirected
                } else if reverse {
                    Orientation::Reversed
                } else {
                    Orientation::Default
                };
                match (edge.source, edge.sink) {
                    (Endpoint::Vertex(source), Endpoint::Vertex(sink)) => KeyedEdge::Paired {
                        source: (NodeIndex(source), None),
                        sink: (NodeIndex(sink), None),
                        key,
                        orientation,
                    },
                    (Endpoint::External(external), Endpoint::Vertex(vertex))
                    | (Endpoint::Vertex(vertex), Endpoint::External(external)) => {
                        let shown = if swapped { partner(external) } else { external };
                        KeyedEdge::Dangling {
                            node: NodeIndex(vertex),
                            hedge_key: Some(shown.index),
                            key,
                            orientation,
                            flow: if edge.source.external().is_some() {
                                Flow::Sink
                            } else {
                                Flow::Source
                            },
                        }
                    }
                    (Endpoint::External(_), Endpoint::External(_)) => {
                        unreachable!("generated edges touch an interaction vertex")
                    }
                }
            })
        };
        let mut candidates = Vec::new();
        for swapped in [false, true] {
            if swapped && !include_left_right_remapping {
                continue;
            }
            let labelling = CanonicalForm::label_edges(
                std::iter::repeat_n(NoData {}, graph.n_vertices()),
                reduced_edges(swapped),
            )
            .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
            candidates.push((labelling, swapped));
        }
        let (reduced, left_right_swapped) = candidates
            .into_iter()
            .min_by(|left, right| (&left.0.form, left.1).cmp(&(&right.0.form, right.1)))
            .ok_or_else(|| GenerationError::GraphConstruction("no canonical graph".to_owned()))?;
        // Reduced colors deliberately forget particle and flow distinctions for numerator
        // grouping. When the reduced group moves a vertex, fix only that residual labelling
        // ambiguity in the physical frame. The choice above compares reduced forms alone, so
        // only the selected exchange needs it.
        let positions = if reduced.moves_nodes {
            Self::anchored_positions(
                graph.vertices(),
                physical_edges(left_right_swapped),
                reduced_edges(left_right_swapped),
                &reduced.form,
            )?
        } else {
            reduced.nodes.into_iter().map(|(_, node)| node.0).collect()
        };

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
        let mut edges = edges
            .into_iter()
            .enumerate()
            .map(|(edge_id, (edge, _, _, data))| {
                // All external edges point into the interaction graph and internal
                // edges point from the lower to the higher canonical vertex.
                let flipped = edge.sink.external().is_some()
                    || matches!(
                        (edge.source, edge.sink),
                        (Endpoint::Vertex(source), Endpoint::Vertex(sink))
                            if positions[source] > positions[sink]
                    );
                // CP acts on neutral adjoint Faddeev-Popov ghosts through color, without
                // exchanging ghost and antighost fields. Their derivative momentum keeps its
                // ghost-number orientation.
                let conjugated = left_right_swapped
                    && !(data.is_ghost() && data.color == 8 && data.charge.is_zero());
                let mut particle = edge.particle;
                if edge.directed && flipped != conjugated {
                    particle = data.antiparticle;
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

    /// The canonical positions of the `reduced` skeleton (canonical form `form`) in an invariant
    /// physical input order. Equivalent labellings are chosen through the node order alone, so
    /// the reduced graph is labelled again with its nodes renumbered into the canonical order of
    /// the `physical` graph on the same vertices (`rules`; particles, orientations and external
    /// states), and its edges keep their order. No reduced color or topology certificate changes,
    /// and no group is expanded. The physical frame depends on more than the reduced input. The
    /// search's choice among equivalent labellings follows the node order through its cell order
    /// and first best leaf, not through the group alone, so neither the first labelling nor its
    /// generators determine the anchored one.
    fn anchored_positions(
        rules: impl ExactSizeIterator<Item = VertexRuleId>,
        physical: impl IntoIterator<Item = KeyedEdge<ParticleId, Option<usize>>>,
        reduced: impl IntoIterator<Item = KeyedEdge<usize, Option<i64>>>,
        form: &CanonicalForm<NoData, usize, Option<i64>>,
    ) -> Result<Vec<usize>, GenerationError> {
        let n_vertices = rules.len();
        let physical = CanonicalForm::label_edges(rules, physical)
            .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?
            .nodes;
        let renumber = |node: NodeIndex| physical[node];
        let anchored = CanonicalForm::label_edges(
            std::iter::repeat_n(NoData {}, n_vertices),
            reduced.into_iter().map(|edge| match edge {
                KeyedEdge::Paired {
                    source,
                    sink,
                    key,
                    orientation,
                } => KeyedEdge::Paired {
                    source: (renumber(source.0), source.1),
                    sink: (renumber(sink.0), sink.1),
                    key,
                    orientation,
                },
                KeyedEdge::Dangling {
                    node,
                    hedge_key,
                    key,
                    orientation,
                    flow,
                } => KeyedEdge::Dangling {
                    node: renumber(node),
                    hedge_key,
                    key,
                    orientation,
                    flow,
                },
            }),
        )
        .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
        debug_assert_eq!(form, &anchored.form);
        Ok(physical
            .into_iter()
            .map(|(_, canonical)| anchored.nodes[canonical].0)
            .collect())
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
        options: &Request,
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

        let blob_range = options.cuts.blobs;
        let spectator_range = options.cuts.spectators;

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
            blob_range.contains(blobs) && spectator_range.contains(spectators)
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
        // Mirror merging compares inventories before finalization, so equal graphs need
        // equally ordered diagram-level cut and threshold keys.
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

    fn passes_sewn_filter(&self, graph: &GeneratedGraph, options: &Request) -> bool {
        if !options.cuts.sewn_tadpoles {
            return true;
        }

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
    use crate::{
        Bridges, CancellationToken, CutAmplitude, Cuts, GraphGroupingOptions, Grouping, Loops,
        Orders, Range, Restrict, Symmetrize, TopologyRules, VertexSelector,
    };
    use DiagramEndpoint::Source;
    use DiagramEndpoint::Target;
    use feynkit_model::Model;
    use idenso::CookMode;
    use idenso::CookSettings;
    use idenso::IndexTooling;
    use idenso::color::ColorSimplifySettings;
    use idenso::representations::ColorAdjoint;
    use idenso::tensor::AlgebraSettings;
    use idenso::tensor::SymbolicTensor;
    use idenso::tensor::composition::TensorCompositionError;
    use idenso::tensor::inference::TensorInferenceError;
    use linnet::half_edge::generate::batch::{Hook, HookMask};
    use serde_json::json;
    use spenso::network::ExecutionResult;
    use spenso::network::Sequential;
    use spenso::network::SmallestDegree;
    use spenso::network::library::symbolic::ExplicitKey;
    use spenso::network::parsing::ParseSettings as TensorParseSettings;
    use spenso::network::parsing::StrictTensorFilter;
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::tensor_symbol;
    use spenso::tensors::parametric::MixedTensor;
    use spenso_hep_lib::FUN_LIB;
    use spenso_hep_lib::HepNet;
    use spenso_hep_lib::hep_lib_atom;
    use std::env;
    use std::process::Command;
    use symbolica::function;
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
        expression.visitor(&mut |term| {
            if matches!(term, AtomView::Fun(function) if function.get_symbol() == head) {
                count += 1;
            }
            true
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

    fn edge(model: &Model, pdg: i64) -> ParticleId {
        particle(model, pdg)
    }

    fn quark_cycle_reflection(model: &Model, reflected: bool) -> GeneratedGraph {
        // The undirected skeleton has a reflection, but the two oriented quark cycles
        // make the physical graph asymmetric. This is the qcd_vac1pi_L5 D2487/D2439 pair.
        let rules = [
            "V_76", "V_76", "V_76", "V_76", "V_75", "V_75", "V_75", "V_36",
        ]
        .map(|name| vertex(model, name));
        let edges = [
            (1, 0, 5),
            (0, 3, 5),
            (0, 6, 21),
            (2, 1, 5),
            (1, 5, 21),
            (3, 2, 5),
            (2, 7, 21),
            (3, 7, 21),
            (4, 5, 3),
            (6, 4, 3),
            (4, 7, 21),
            (5, 6, 3),
        ]
        .map(|(mut source, mut sink, pdg)| {
            let directed = pdg != 21;
            if directed && reflected {
                std::mem::swap(&mut source, &mut sink);
            }
            GeneratedEdge {
                source: Endpoint::Vertex(source),
                sink: Endpoint::Vertex(sink),
                particle: particle(model, pdg),
                directed,
            }
        });
        GeneratedGraph::new(rules, edges)
    }

    #[test]
    fn canonical_physical_frame_resolves_reduced_skeleton_reflections() {
        let model = Arc::new(Model::qcd());
        let process = Process::new([] as [i64; 0], [] as [i64; 0]);
        let request = process.clone().amplitude(Loops::total(5..=5));
        let resolved =
            ResolvedProcess::new(&model, &process, GenerationType::Amplitude, &request).unwrap();
        let generator = Generator::new(Arc::clone(&model));
        let graphs = [false, true].map(|reflected| quark_cycle_reflection(&model, reflected));
        assert_eq!(
            graphs[0].to_hedge_graph().canonical_form(),
            graphs[1].to_hedge_graph().canonical_form(),
            "the input pair is the same directed, rule-coloured physical graph"
        );
        let mut variants = graphs.to_vec();
        let permutation = [3, 5, 0, 7, 1, 6, 2, 4];
        let mut rules = vec![vertex(&model, "V_36"); permutation.len()];
        for (old, rule) in graphs[0].vertices().enumerate() {
            rules[permutation[old]] = rule;
        }
        let edges = graphs[0]
            .edges()
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .map(|mut edge| {
                let relabel = |endpoint| match endpoint {
                    Endpoint::Vertex(vertex) => Endpoint::Vertex(permutation[vertex]),
                    external => external,
                };
                edge.source = relabel(edge.source);
                edge.sink = relabel(edge.sink);
                if edge.directed {
                    std::mem::swap(&mut edge.source, &mut edge.sink);
                    edge.particle = model.particle_by_id(edge.particle).unwrap().antiparticle;
                }
                edge
            })
            .collect::<Vec<_>>();
        variants.push(GeneratedGraph::new(rules, edges));
        // A sibling decoration over the first graph's incidence, every quark replaced by its
        // antiquark: physically the reflected graph, with the same reduced canonizer input
        // when only masses are compared. Its canonical frame places the vertices differently,
        // so no labelling may be shared through that input or the incidence.
        let sibling = graphs[0].reframed(
            &(0..graphs[0].n_vertices()).collect::<Vec<_>>(),
            graphs[0].edges().enumerate().map(|(edge, data)| {
                let particle = if data.directed {
                    model.particle_by_id(data.particle).unwrap().antiparticle
                } else {
                    data.particle
                };
                (edge, false, particle)
            }),
            |external| external,
        );
        for grouping in [
            Grouping::None,
            Grouping::Identical(GraphGroupingOptions::default()),
        ] {
            let positions = |graph: &GeneratedGraph| {
                let canonical = resolved
                    .canonicalize_numerator_graph(graph, &model, &grouping)
                    .unwrap();
                graph.frame_map(&canonical).vertices
            };
            assert_ne!(positions(&graphs[0]), positions(&sibling));
        }
        variants.push(sibling);
        let numerators = NumeratorOptions::default();
        for grouping in [
            Grouping::None,
            Grouping::Identical(GraphGroupingOptions::default()),
            Grouping::Identical(GraphGroupingOptions {
                differentiate_particle_masses_only: false,
                ..GraphGroupingOptions::default()
            }),
        ] {
            let diagrams = variants
                .iter()
                .map(|graph| {
                    let input = generator
                        .diagram_input(&resolved, graph.clone(), 1, 1, &request, &grouping)
                        .unwrap()
                        .unwrap();
                    generator
                        .diagram_pair(
                            Some("D".into()),
                            GenerationType::Amplitude,
                            input,
                            (&numerators, &numerators),
                        )
                        .unwrap()
                        .into_representative()
                        .unwrap()
                })
                .collect::<Vec<_>>();
            for other in &diagrams[1..] {
                assert_eq!(
                    diagrams[0].canonical_key().unwrap(),
                    other.canonical_key().unwrap(),
                    "physical slot assignment must respect the typed graph isomorphism"
                );
                assert_eq!(
                    serde_json::to_vec(&diagrams[0]).unwrap(),
                    serde_json::to_vec(other).unwrap(),
                    "isomorphic representatives need the same numerator and routing frame"
                );
            }
        }
    }

    #[test]
    fn canonical_physical_frame_accepts_reversed_external_charge_representations() {
        let model = Arc::new(Model::qcd());
        let process = Process::new([2_i64, -2], [] as [i64; 0]);
        let request = process.clone().amplitude(Loops::total(5..=5));
        let resolved =
            ResolvedProcess::new(&model, &process, GenerationType::Amplitude, &request).unwrap();
        let generator = Generator::new(Arc::clone(&model));
        let base = quark_cycle_reflection(&model, false);
        let mut rules = base.vertices().collect::<Vec<_>>();
        rules.extend([vertex(&model, "V_36"), vertex(&model, "V_135")]);
        let mut edges = base.edges().collect::<Vec<_>>();
        // Attach an external quark current to the fixed gluon branch of the reflection.
        edges[10].sink = Endpoint::Vertex(8);
        edges.extend([(8, 7), (8, 9)].map(|(source, sink)| GeneratedEdge {
            source: Endpoint::Vertex(source),
            sink: Endpoint::Vertex(sink),
            particle: particle(&model, 21),
            directed: false,
        }));
        for (index, pdg) in [2, -2].into_iter().enumerate() {
            let external = Endpoint::External(ExternalNode {
                index,
                state: ExternalState::Incoming,
                particle: particle(&model, pdg),
            });
            let (source, sink) = if pdg > 0 {
                (external, Endpoint::Vertex(9))
            } else {
                (Endpoint::Vertex(9), external)
            };
            edges.push(GeneratedEdge {
                source,
                sink,
                particle: particle(&model, 2),
                directed: true,
            });
        }
        let reversed = edges.iter().rev().copied().map(|mut edge| {
            if edge.directed {
                std::mem::swap(&mut edge.source, &mut edge.sink);
                edge.particle = model.particle_by_id(edge.particle).unwrap().antiparticle;
            }
            edge
        });
        let graphs = [
            GeneratedGraph::new(rules.iter().copied(), edges.iter().copied()),
            GeneratedGraph::new(rules.iter().copied(), reversed),
        ];
        let numerators = NumeratorOptions::default();
        let diagrams = graphs.map(|graph| {
            let input = generator
                .diagram_input(&resolved, graph, 1, 1, &request, &Grouping::None)
                .unwrap()
                .unwrap();
            generator
                .diagram_pair(
                    Some("D".into()),
                    GenerationType::Amplitude,
                    input,
                    (&numerators, &numerators),
                )
                .unwrap()
                .into_representative()
                .unwrap()
        });
        assert_eq!(
            serde_json::to_vec(&diagrams[0]).unwrap(),
            serde_json::to_vec(&diagrams[1]).unwrap()
        );
    }

    #[test]
    fn canonical_physical_frame_preserves_pure_edge_kernels() {
        let model = Model::phi4();
        let phi = model.particle_id("phi").unwrap();
        let rule = model.vertex_rule_id_at(0).unwrap();
        let process = Process::new([] as [i64; 0], [] as [i64; 0]);
        let request = process.clone().amplitude(Loops::total(2..=2));
        let resolved =
            ResolvedProcess::new(&model, &process, GenerationType::Amplitude, &request).unwrap();
        let graph = GeneratedGraph::new(
            [rule],
            [GeneratedEdge {
                source: Endpoint::Vertex(0),
                sink: Endpoint::Vertex(0),
                particle: phi,
                directed: false,
            }; 2],
        );
        let hedge_graph = graph.to_hedge_graph();
        let canonical = hedge_graph.canonize();
        assert_eq!(canonical.automorphism_count(), Some(8));
        assert!(
            canonical
                .generators()
                .all(|g| g.node(NodeIndex(0)) == NodeIndex(0))
        );
        assert_eq!(
            resolved
                .canonicalize_numerator_graph(&graph, &model, &Grouping::None)
                .unwrap(),
            graph
        );
    }

    #[test]
    fn canonical_physical_frame_preserves_external_labels_and_flavour_reflections() {
        let model = Arc::new(Model::qcd());
        let process = Process::new([21_i64], [21_i64]);
        let request = process.clone().amplitude(Loops::total(4..=4));
        let resolved =
            ResolvedProcess::new(&model, &process, GenerationType::Amplitude, &request).unwrap();
        let generator = Generator::new(Arc::clone(&model));
        // qcd_gse_L4: two equal-mass quark triangles and fixed external gluons.
        let rules = [
            "V_135", "V_75", "V_75", "V_135", "V_75", "V_135", "V_36", "V_36",
        ]
        .map(|name| vertex(&model, name));
        let mut edges = [
            (0, 6, ExternalState::Incoming),
            (1, 7, ExternalState::Outgoing),
        ]
        .map(|(index, vertex, state)| GeneratedEdge {
            source: Endpoint::External(ExternalNode {
                index,
                state,
                particle: particle(&model, 21),
            }),
            sink: Endpoint::Vertex(vertex),
            particle: particle(&model, 21),
            directed: false,
        })
        .to_vec();
        edges.extend(
            [
                (0, 1, 21),
                (3, 0, 2),
                (0, 5, 2),
                (1, 2, 3),
                (4, 1, 3),
                (2, 4, 3),
                (2, 7, 21),
                (5, 3, 2),
                (3, 7, 21),
                (4, 6, 21),
                (5, 6, 21),
            ]
            .map(|(source, sink, pdg)| GeneratedEdge {
                source: Endpoint::Vertex(source),
                sink: Endpoint::Vertex(sink),
                particle: particle(&model, pdg),
                directed: pdg != 21,
            }),
        );
        let permutation = [1, 0, 3, 2, 5, 4, 6, 7];
        let mut reflected_rules = rules;
        for (vertex, rule) in rules.into_iter().enumerate() {
            reflected_rules[permutation[vertex]] = rule;
        }
        let reflected_edges = edges.iter().copied().map(|mut edge| {
            let reflect = |endpoint| match endpoint {
                Endpoint::Vertex(vertex) => Endpoint::Vertex(permutation[vertex]),
                external => external,
            };
            edge.source = reflect(edge.source);
            edge.sink = reflect(edge.sink);
            edge
        });
        let graphs = [
            GeneratedGraph::new(rules, edges.iter().copied()),
            GeneratedGraph::new(reflected_rules, reflected_edges),
        ];
        let numerators = NumeratorOptions::default();
        for exact_species in [false, true] {
            let grouping = Grouping::Identical(GraphGroupingOptions {
                differentiate_particle_masses_only: !exact_species,
                ..GraphGroupingOptions::default()
            });
            let diagrams = graphs.each_ref().map(|graph| {
                let input = generator
                    .diagram_input(&resolved, graph.clone(), 1, 1, &request, &grouping)
                    .unwrap()
                    .unwrap();
                generator
                    .diagram_pair(
                        Some("D".into()),
                        GenerationType::Amplitude,
                        input,
                        (&numerators, &numerators),
                    )
                    .unwrap()
                    .into_representative()
                    .unwrap()
            });
            assert_eq!(
                serde_json::to_vec(&diagrams[0]).unwrap(),
                serde_json::to_vec(&diagrams[1]).unwrap(),
                "external labels, local slots and physical frames must remain equivariant"
            );
        }
    }

    #[test]
    fn canonical_physical_frame_preserves_group_members_and_ratios() {
        let model = Arc::new(Model::qed());
        let case = cases::Case::named("qed_eeaa_scalar_L1");
        let run = Run::default().threads(1);
        for masses_only in [true, false] {
            for group in [
                Grouping::Identical,
                Grouping::UpToSign,
                Grouping::UpToScalar,
            ] {
                let grouping = group(GraphGroupingOptions {
                    differentiate_particle_masses_only: masses_only,
                    ..GraphGroupingOptions::default()
                });
                let mut reference = None;
                for classes in [crate::Classes::None, crate::Classes::Auto] {
                    let built = case
                        .request(&model)
                        .with_classes(classes)
                        .generate(Arc::clone(&model), &run)
                        .unwrap()
                        .build(
                            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                            &run,
                        )
                        .unwrap();
                    let grouped = built.group(&grouping, &run).unwrap();
                    let groups = grouped
                        .groups()
                        .iter()
                        .map(|group| {
                            let members = group
                                .members
                                .iter()
                                .map(|member| {
                                    (
                                        member.source_identity.as_ref().unwrap().id,
                                        (
                                            member.ratio.to_canonical_string(),
                                            member.overall_factor.to_canonical_string(),
                                        ),
                                    )
                                })
                                .collect::<BTreeMap<_, _>>();
                            assert_eq!(members.len(), group.members.len());
                            (grouped.get(group.master).unwrap().id().unwrap(), members)
                        })
                        .collect::<BTreeMap<_, _>>();
                    assert_eq!(groups.len(), grouped.groups().len());
                    if let Some(reference) = &reference {
                        assert_eq!(&groups, reference, "{grouping:?}: exact members and ratios");
                    } else {
                        reference = Some(groups);
                    }
                }
            }
        }
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub(super) enum ColoredNode {
        External(ExternalNode),
        Interaction(VertexRuleId),
    }

    /// A generated-graph fixture written with degree-one external vertices; [`Self::build`]
    /// turns them into legs and keeps the edge order.
    #[derive(Debug, Clone, Default)]
    pub(super) struct Graph {
        nodes: Vec<ColoredNode>,
        edges: Vec<(usize, usize, bool, ParticleId)>,
    }

    impl Graph {
        pub(super) fn new() -> Self {
            Self::default()
        }

        pub(super) fn add_node(&mut self, node: ColoredNode) -> usize {
            self.nodes.push(node);
            self.nodes.len() - 1
        }

        /// Undirected edges are stored from the lower to the higher node, as Graphica did.
        pub(super) fn add_edge(
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

        pub(super) fn build(&self) -> GeneratedGraph {
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
        let automatic = process
            .clone()
            .amplitude(Loops::total(0..=1))
            .with_topology(
                TopologyRules::new()
                    .with_self_loops(false)
                    .with_max_vertices(4)
                    .with_bridges(Bridges::OnePIInLoops),
            );
        assert_eq!(automatic.topology().bridges, Bridges::OnePIInLoops);

        let all = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .amplitude(Loops::total(0..=1))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(4),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        let filtered = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .amplitude(Loops::total(0..=1))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(4)
                        .with_bridges(Bridges::OnePIInLoops),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        let strict = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .amplitude(Loops::total(0..=1))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(4)
                        .with_bridges(Bridges::AtMost(0)),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        let ids = |result: &BuiltDiagramSet, loops| {
            result
                .diagrams()
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

        let loops_only = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .amplitude(Loops::total(1..=1))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(4)
                        .with_bridges(Bridges::OnePIInLoops),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
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
        let generated = {
            let run = Run::default();
            let request = (process)
                .clone()
                .amplitude(Loops::total(0..=1))
                .with_topology(TopologyRules::new().with_self_loops(true));
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        assert!(!generated.diagrams().is_empty());
        assert!(generated.report().completed);
        assert!(
            generated
                .diagrams()
                .iter()
                .all(|diagram| !diagram.numerator().unwrap().is_zero())
        );
    }

    #[test]
    fn finalizes_external_fermion_projectors_unless_explicitly_overridden() {
        let generator = Generator::new(fermion_model());
        let process = Process::new([1_i64, -1], [1_i64, -1]);

        let generated = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .amplitude(Loops::total(0..=0))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(2),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        assert!(!generated.diagrams().is_empty());
        for diagram in generated.diagrams() {
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

        let explicit = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .amplitude(Loops::total(0..=0))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(2),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default()
                            .with_naming(Naming::prefix("D"))
                            .with_projector(Atom::one()),
                        &run,
                    )
                })
        }
        .unwrap();
        assert!(
            explicit
                .diagrams()
                .iter()
                .all(|diagram| diagram.projector().expression() == &Atom::one())
        );
    }

    #[test]
    fn compton_cuts_preserve_the_requested_particle_charge() {
        let generator = Generator::new(fermion_model());
        for pdg in [1_i64, -1] {
            let process = Process::new([pdg, 22], [pdg, 22]);
            let generated = {
                let run = Run::default().threads(1);
                let request = (process)
                    .clone()
                    .cross_section(Loops::forward(1..=1))
                    .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                    .with_topology(
                        TopologyRules::new()
                            .with_self_loops(false)
                            .with_max_vertices(4),
                    );
                request
                    .generate(generator.model.clone(), &run)
                    .and_then(|prepared| {
                        prepared.build(
                            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                            &run,
                        )
                    })
            }
            .unwrap();
            assert!(!generated.diagrams().is_empty());
            for diagram in generated.diagrams() {
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
            let generated = {
                let run = Run::default().threads(1);
                let request = (process)
                    .clone()
                    .cross_section(Loops::forward(1..=1))
                    .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                    .with_topology(
                        TopologyRules::new()
                            .with_self_loops(false)
                            .with_max_vertices(4),
                    );
                request
                    .generate(generator.model.clone(), &run)
                    .and_then(|prepared| {
                        prepared.build(
                            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                            &run,
                        )
                    })
            }
            .unwrap();
            assert!(!generated.diagrams().is_empty());
            for diagram in generated.diagrams() {
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
        let generated = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .cross_section(Loops::forward(0..=1))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(4),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        assert!(!generated.diagrams().is_empty());

        for diagram in generated.diagrams() {
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
        for self_energy_on_source_side in [None, Some(true), Some(false)] {
            let model = fermion_model();
            let resolved = ResolvedProcess::new(
                &model,
                &Process::new([22_i64, 22], [22_i64, 22]),
                GenerationType::CrossSection,
                &(Process::new([22_i64, 22], [22_i64, 22]))
                    .clone()
                    .cross_section(Loops::forward(0..=0))
                    .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                    .with_topology(TopologyRules::new().with_self_loops(false)),
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
            for (source_side, (a, b)) in [(true, (left_0, left_1)), (false, (right_0, right_1))] {
                if self_energy_on_source_side == Some(source_side) {
                    // The seeded attachments enclose this pocket on one side.
                    // Cut enumeration must absorb it before testing connectivity.
                    let m0 = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
                    let m1 = graph.add_node(ColoredNode::Interaction(vertex(&model, "V")));
                    graph.add_edge(a, m0, false, photon).unwrap();
                    graph.add_edge(m0, m1, false, photon).unwrap();
                    graph.add_edge(m0, m1, false, photon).unwrap();
                    graph.add_edge(m1, b, false, photon).unwrap();
                } else {
                    graph.add_edge(a, b, false, photon).unwrap();
                }
            }

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
                .cut_partitions(
                    &graph.build(),
                    &model,
                    &(Process::new([22_i64, 22], [22_i64, 22]))
                        .clone()
                        .cross_section(Loops::forward(0..=0))
                        .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                        .with_topology(TopologyRules::new().with_self_loops(false)),
                )
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

                    numerator: SymbolicTensor::infer(Atom::one()).unwrap(),
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
            let process = Process::new(incoming, outgoing);
            for symmetrize in [false, true] {
                let run = Run::default().threads(1);
                let request = process
                    .clone()
                    .cross_section(Loops::forward(loops..=loops))
                    .with_topology(TopologyRules::new().with_self_loops(false))
                    .with_restrict(Restrict::new().only_particles(selection.iter().copied()))
                    .with_symmetrize(Symmetrize {
                        initial: symmetrize,
                        ..Symmetrize::default()
                    })
                    .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0));
                let result = request
                    .generate(model.clone().into(), &run)
                    .unwrap()
                    .build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                    .unwrap();
                assert!(!result.diagrams().is_empty(), "{incoming:?} > {outgoing:?}");
                let several_cuts = result
                    .diagrams()
                    .iter()
                    .any(|diagram| diagram.cuts().len() > 1);
                assert_eq!(several_cuts, loops > 1, "{incoming:?} > {outgoing:?}");
                for diagram in result.diagrams() {
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
        let vectors = {
            let run = Run::default().threads(1);
            let request = (photon_process())
                .clone()
                .amplitude(Loops::total(1..=1))
                .with_topology(TopologyRules::new().with_self_loops(false));
            request
                .generate((fermion_model()).into(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        assert!(!vectors.diagrams().is_empty());
        for diagram in vectors.diagrams() {
            for (head, count) in [
                (feynkit_graph::symbols::epsilon(), 1),
                (feynkit_graph::symbols::epsilonbar(), 1),
                (feynkit_graph::symbols::hedge_index(), 2),
            ] {
                assert_eq!(
                    function_count(&NumeratorInterface::export(diagram.projector()).0, head),
                    count
                );
            }
        }

        let scalars = {
            let run = Run::default().threads(1);
            let request = (Process::new(["phi"], ["phi"]))
                .clone()
                .amplitude(Loops::total(1..=1))
                .with_topology(TopologyRules::new().with_self_loops(false));
            request
                .generate((scalar_model()).into(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        assert!(!scalars.diagrams().is_empty());
        assert!(
            scalars
                .diagrams()
                .iter()
                .all(|diagram| diagram.projector().expression() == &Atom::one())
        );
    }

    #[test]
    fn forced_cuts_retain_only_the_requested_typed_edge_set() {
        let generator = Generator::new(scalar_model());
        let process = scalar_cross_section_process();

        let generated = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .cross_section(Loops::forward(1..=1))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(4),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        let requested: BTreeSet<EdgeId> = generated
            .diagrams()
            .iter()
            .flat_map(|diagram| diagram.cuts())
            .next()
            .unwrap()
            .cut
            .iter()
            .map(|half_edge| half_edge.edge)
            .collect();

        let forced = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .cross_section(Loops::forward(1..=1))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(4),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
                .and_then(|built| {
                    built.select(
                        &Selection::new().with_forced_cuts([requested.clone()]),
                        &run,
                    )
                })
        }
        .unwrap();

        assert!(!forced.diagrams().is_empty());
        assert!(forced.diagrams().iter().all(|diagram| {
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
        let result = {
            let run = Run::default().threads(1);
            let request = (scalar_cross_section_process())
                .clone()
                .cross_section(Loops::forward(1..=1))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(4),
                );
            request
                .generate((scalar_model()).into(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
                .and_then(|built| built.select(&Selection::new().with_forced_cuts([unknown]), &run))
        };

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

        let generated = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .cross_section(Loops::forward(1..=1))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(4),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        let cut_sets = generated
            .diagrams()
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

        let grouped = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .cross_section(Loops::forward(1..=1))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(4),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
                .and_then(|built| {
                    built.select(
                        &Selection::new().with_forced_cuts([requested.clone()]),
                        &run,
                    )
                })
                .and_then(|built| {
                    built.group(&Grouping::Identical(GraphGroupingOptions::default()), &run)
                })
        }
        .unwrap();

        grouped.validate_groups().unwrap();
        assert!(!grouped.groups().is_empty());
        assert!(grouped.diagrams().iter().all(|diagram| {
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
            {
                let run = Run::default();
                let request = (process)
                    .clone()
                    .amplitude(Loops::total(0..=0))
                    .with_topology(TopologyRules::new().with_self_loops(false));
                request
                    .generate((target.clone()).into(), &run)
                    .and_then(|prepared| {
                        prepared.build(
                            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                            &run,
                        )
                    })
            },
            Err(GenerationError::Selector(SelectorError::ModelMismatch {
                kind: "particle",
                ..
            }))
        ));

        let foreign_vertex =
            VertexSelector::by_id(&source, source.vertex_rule_id("V1").unwrap()).unwrap();
        for filter in [
            Restrict::new().exclude_particles([foreign_particle]),
            Restrict::new().only_vertices([foreign_vertex]),
        ] {
            let request = Process::new(["phi"], ["phi", "phi"])
                .amplitude(Loops::total(0..=0))
                .with_restrict(filter);
            assert!(matches!(
                request.generate(target.clone().into(), &Run::default()),
                Err(GenerationError::Selector(
                    SelectorError::ModelMismatch { .. }
                ))
            ));
        }
    }

    #[test]
    fn typed_search_hooks_separate_complete_dressing_from_selection() {
        let generator = Generator::new(scalar_model());
        let skeletons = Arc::new(AtomicUsize::new(0));
        let dressings = Arc::new(AtomicUsize::new(0));
        let selections = Arc::new(AtomicUsize::new(0));
        let request = Process::new(["phi", "phi"], ["phi", "phi"])
            .amplitude(Loops::total(0..=0))
            .with_topology(TopologyRules::new().with_self_loops(false))
            .prune_skeletons({
                let count = Arc::clone(&skeletons);
                move |view| {
                    count.fetch_add(usize::from(view.is_topology_complete()), Ordering::Relaxed);
                    true
                }
            })
            .prune_dressings({
                let count = Arc::clone(&dressings);
                move |view| {
                    count.fetch_add(usize::from(view.is_complete()), Ordering::Relaxed);
                    true
                }
            })
            .select({
                let count = Arc::clone(&selections);
                move |view| {
                    assert!(view.is_complete());
                    count.fetch_add(1, Ordering::Relaxed);
                    false
                }
            });
        let run = Run::default().threads(1);
        let result = request
            .generate(generator.model.clone(), &run)
            .unwrap()
            .build(
                &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                &run,
            )
            .unwrap();
        assert!(result.report().completed);
        assert!(result.diagrams().is_empty());
        assert!(skeletons.load(Ordering::Relaxed) > 0);
        assert!(dressings.load(Ordering::Relaxed) > 0);
        assert_eq!(
            dressings.load(Ordering::Relaxed),
            selections.load(Ordering::Relaxed)
        );
    }

    #[test]
    fn deferred_typed_hooks_preserve_outputs_and_cancel_without_worker_replies() {
        let generator = Generator::new(scalar_model());
        let process = Process::new(["phi", "phi"], ["phi", "phi", "phi"]);
        let request = process
            .clone()
            .amplitude(Loops::total(0..=0))
            .with_topology(TopologyRules::new().with_self_loops(false));
        let expected = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .amplitude(Loops::total(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false));
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        let json = |result: &BuiltDiagramSet| {
            result
                .diagrams()
                .iter()
                .map(|diagram| diagram.to_json().unwrap())
                .collect::<Vec<_>>()
        };
        let mask = HookMask {
            skeletons: true,
            dressings: true,
            select: true,
        };
        for threads in [1, 4] {
            let observed = Arc::new(Mutex::new([0usize; 3]));
            let callbacks = Arc::clone(&observed);
            let run = Run::default().threads(threads);
            let deferred = request.clone().deferred_hooks(mask, move |requests| {
                let mut observed = callbacks.lock().unwrap();
                ControlFlow::Continue(
                    requests
                        .into_iter()
                        .map(|request| {
                            match request.hook() {
                                Hook::Skeleton => {
                                    assert!(request.skeleton().is_some());
                                    assert!(request.dressing().is_none());
                                    observed[0] += 1;
                                }
                                Hook::Dressing => {
                                    let view = request.dressing().unwrap();
                                    observed[1] += usize::from(view.is_complete());
                                }
                                Hook::Select => {
                                    assert!(request.diagram().unwrap().is_complete());
                                    observed[2] += 1;
                                }
                            }
                            request.reply(true)
                        })
                        .collect(),
                )
            });
            let actual = deferred
                .generate(generator.model.clone(), &run)
                .unwrap()
                .build(
                    &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                    &run,
                )
                .unwrap();
            assert!(actual.report().completed);
            assert_eq!(json(&actual), json(&expected), "{threads} threads");
            assert_eq!(
                actual.report().topology_count,
                expected.report().topology_count
            );
            let observed = observed.lock().unwrap();
            assert!(observed.iter().all(|&count| count > 0));
            assert_eq!(observed[1], observed[2]);
        }
        let cancelled = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .amplitude(Loops::total(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false))
                .deferred_hooks(mask, |_| ControlFlow::Break(()));
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        assert!(!cancelled.report().completed);
        assert!(cancelled.diagrams().is_empty());
    }

    #[test]
    fn empty_generation_returns_a_completed_result() {
        let generator = Generator::new(scalar_model());
        let process = Process::new(["phi"], ["phi", "phi"]);
        let request = process
            .amplitude(Loops::total(0..=0))
            .with_topology(TopologyRules::new().with_self_loops(false));
        for request in [
            request.clone().with_topology(
                TopologyRules::new()
                    .with_self_loops(false)
                    .with_max_vertices(0),
            ),
            request.with_restrict(Restrict::new().exclude_particles(["phi"])),
        ] {
            let run = Run::default();
            let result = request
                .generate(generator.model.clone(), &run)
                .unwrap()
                .build(
                    &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                    &run,
                )
                .unwrap();
            assert!(result.report().completed);
            assert!(result.diagrams().is_empty());
            assert!(result.groups().is_empty());
            assert_eq!(result.report().topology_count, 0usize.into());
            assert_eq!(result.report().interaction_assignment_count, 0usize.into());
            assert_eq!(result.report().retained_count, 0usize.into());
        }
    }

    #[test]
    fn cancelled_generation_returns_an_incomplete_result() {
        let cancellation = CancellationToken::new();
        cancellation.cancel();
        let generated = {
            let run = Run::default().cancellation_token(cancellation);
            let request = (Process::new(["phi"], ["phi", "phi"]))
                .clone()
                .amplitude(Loops::total(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false));
            request
                .generate((scalar_model()).into(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();

        assert!(!generated.report().completed);
        assert!(generated.diagrams().is_empty());
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
        counterterm["name"] = json!("V1_CT");
        rules.push(counterterm);
        let model: Model = serde_json::from_value(definition).unwrap();
        let generated = {
            let run = Run::default().threads(1);
            let request = (Process::new(["phi"], ["phi"]))
                .clone()
                .amplitude(Loops::total(1..=1))
                .with_topology(TopologyRules::new().with_self_loops(false));
            request
                .generate((model.clone()).into(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        // The bubble with either rule at either vertex: the legs distinguish the vertices, and
        // the two internal lines are interchangeable.
        let mut rules = generated
            .diagrams()
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

        let ids = |threads| {
            {
                let run = Run::default().threads(threads);
                let request = (process)
                    .clone()
                    .amplitude(Loops::total(0..=0))
                    .with_topology(TopologyRules::new().with_self_loops(false))
                    .with_symmetrize(Symmetrize {
                        initial: true,
                        final_state: true,
                        ..Symmetrize::default()
                    });
                request
                    .generate(generator.model.clone(), &run)
                    .and_then(|prepared| {
                        prepared.build(
                            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                            &run,
                        )
                    })
            }
            .unwrap()
            .diagrams()
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
        let run = Run::default();
        let prepared = process
            .amplitude(Loops::total(0..=0))
            .with_topology(TopologyRules::new().with_self_loops(false))
            .generate(model.into(), &run)
            .unwrap();
        let minimal = NumeratorOptions::default();
        let all = prepared.build(&minimal, &run).unwrap();
        assert_eq!(all.diagrams().len(), 2);
        assert!(
            all.diagrams()
                .iter()
                .all(|diagram| diagram.name().is_none())
        );
        assert!(
            all.groups()
                .iter()
                .flat_map(|group| &group.members)
                .all(|member| member.source_identity.is_none())
        );
        let chosen = all.diagrams()[0].id().unwrap();
        let selected = all
            .select(&Selection::new().select_ids([chosen]), &run)
            .unwrap();
        assert_eq!(selected.diagrams().len(), 1);
        assert_eq!(selected.diagrams()[0].id().unwrap(), chosen);
        assert!(selected.diagrams()[0].name().is_none());
        let vetoed = all
            .select(&Selection::new().veto_ids([chosen]), &run)
            .unwrap();
        assert_eq!(vetoed.diagrams().len(), 1);
        assert_ne!(vetoed.diagrams()[0].id().unwrap(), chosen);
        for empty in [
            Selection::new().select_ids([]),
            Selection::new().select_names(Vec::<String>::new()),
        ] {
            assert!(all.select(&empty, &run).unwrap().diagrams().is_empty());
        }
        let named = minimal.clone().with_naming(Naming::Named {
            prefix: "test_".into(),
        });
        let selected = prepared
            .build(&named, &run)
            .unwrap()
            .select(&Selection::new().select_names(["test_0".to_owned()]), &run)
            .unwrap();
        assert_eq!(selected.diagrams().len(), 1);
        assert_eq!(selected.diagrams()[0].name(), Some("test_0"));
        assert_eq!(selected.diagrams()[0].id().unwrap(), chosen);
        // Both requests match: ID routing takes precedence over the invalid name map.
        let routed = named
            .with_loop_momentum_basis(chosen, [])
            .with_named_loop_momentum_basis("test_0", [EdgeId(usize::MAX)]);
        assert_eq!(prepared.build(&routed, &run).unwrap().diagrams().len(), 2);
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
            definition[field] = json!([]);
        }
        let template = definition["particles"][0].clone();
        definition["particles"] = json!([]);
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
        let process = Process::new(["g"], ["g"]);

        let serial = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .amplitude(Loops::total(2..=2))
                .with_restrict(Restrict::new().exclude_particles(["u", "d", "c", "s", "t"]))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_bridges(Bridges::AtMost(0)),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        let progress = Arc::new(Mutex::new(Vec::new()));
        let progress_sink = Arc::clone(&progress);
        let parallel = {
            let run = Run::default().threads(4).progress(move |update| {
                if update.stage == Stage::Numerators {
                    progress_sink.lock().unwrap().push(update.completed);
                }
                crate::GenerationControl::Continue
            });
            let request = (process)
                .clone()
                .amplitude(Loops::total(2..=2))
                .with_restrict(Restrict::new().exclude_particles(["u", "d", "c", "s", "t"]))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_bridges(Bridges::AtMost(0)),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        assert!(!serial.diagrams().is_empty());
        assert_eq!(serial.report(), parallel.report());
        assert_eq!(
            *progress.lock().unwrap(),
            (0..=parallel.diagrams().len())
                .map(BigUint::from)
                .collect::<Vec<_>>()
        );
        for (serial, parallel) in serial.diagrams().iter().zip(parallel.diagrams()) {
            assert_eq!(serial.to_json().unwrap(), parallel.to_json().unwrap());
        }
        assert_eq!(serial.groups().len(), parallel.groups().len());
        for (serial, parallel) in serial.groups().iter().zip(parallel.groups()) {
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
    fn batched_preparation_report_matches_native_enumeration() {
        let model = Arc::new(scalar_model());
        let request = Process::new(["phi"], ["phi"])
            .amplitude(Loops::total(4..=4))
            .with_topology(
                TopologyRules::new()
                    .with_max_vertices(8)
                    .with_self_loops(false),
            );
        let search = Generator::new(Arc::clone(&model))
            .search(&request, |_, _| {}, |_, _| {})
            .unwrap();
        let assignments = search
            .plan
            .as_ref()
            .unwrap()
            .diagrams(&EngineRun::default())
            .unwrap()
            .count();
        let mut outputs = Vec::new();
        for threads in [1, 4] {
            let progress = Arc::new(Mutex::new(Vec::new()));
            let captured = Arc::clone(&progress);
            let run = Run::default().threads(threads).progress(move |update| {
                if update.stage == Stage::Prepare {
                    captured.lock().unwrap().push(update.completed);
                }
                crate::GenerationControl::Continue
            });
            let prepared = request.generate(Arc::clone(&model), &run).unwrap();
            assert!(prepared.report().completed);
            // Four workers must prepare several bounded batches, rather than only the first
            // batch, and keep the enumeration order.
            assert!(prepared.len() > 32);
            assert_eq!(
                prepared.report().interaction_assignment_count,
                BigUint::from(assignments)
            );
            assert_eq!(
                *progress.lock().unwrap(),
                (1..=assignments).map(BigUint::from).collect::<Vec<_>>()
            );
            outputs.push(prepared.build(&NumeratorOptions::default(), &run).unwrap());
        }
        assert_eq!(outputs[0].report(), outputs[1].report());
        for (serial, parallel) in outputs[0].iter().zip(outputs[1].iter()) {
            assert_eq!(serial.to_json().unwrap(), parallel.to_json().unwrap());
        }
    }

    #[test]
    fn zero_color_filter_is_opt_in_for_generated_qcd_tadpoles() {
        let model = Arc::new(Model::qcd());
        let process = Process::new([] as [&str; 0], ["g"]);
        let request = process.amplitude(Loops::total(1..=1)).with_topology(
            TopologyRules::new()
                .with_max_vertices(1)
                .with_self_loops(true)
                .with_zero_flow_edges(true),
        );
        let run = Run::default().threads(1);
        let numerators = NumeratorOptions::default().with_naming(Naming::prefix("D"));
        let prepared = request.generate(model, &run).unwrap();
        let unfiltered = prepared.build(&numerators, &run).unwrap();
        assert!(unfiltered.report().completed);
        assert!(!unfiltered.diagrams().is_empty());
        assert_eq!(unfiltered.report().zero_numerator_count, 0usize.into());
        assert!(
            unfiltered
                .diagrams()
                .iter()
                .any(|diagram| !diagram.numerator().unwrap().is_zero())
        );

        let updates = Arc::new(Mutex::new(Vec::new()));
        let captured = updates.clone();
        let filtered_run = run.clone().progress(move |update| {
            if update.stage == Stage::ZeroColor {
                captured.lock().unwrap().push(update);
            }
            crate::GenerationControl::Continue
        });
        let filtered = prepared
            .build(
                &numerators.clone().with_zero_color_filter(true),
                &filtered_run,
            )
            .unwrap();
        assert!(filtered.report().completed);
        assert!(filtered.diagrams().is_empty());
        assert_eq!(
            filtered.report().zero_numerator_count,
            unfiltered.diagrams().len().into()
        );
        filtered.validate_groups().unwrap();
        let updates = updates.lock().unwrap();
        assert_eq!(
            updates
                .iter()
                .map(|update| (
                    update.completed.clone(),
                    update.zero_numerator_count.clone()
                ))
                .collect::<Vec<_>>(),
            (0..=unfiltered.diagrams().len())
                .map(|count| (BigUint::from(count), BigUint::from(count)))
                .collect::<Vec<_>>()
        );
        assert!(
            updates
                .iter()
                .all(|update| update.total == Some(unfiltered.diagrams().len().into()))
        );

        let amplitude = unfiltered
            .clone()
            .into_amplitude(Default::default())
            .unwrap();
        assert_eq!(amplitude.diagrams().len(), unfiltered.diagrams().len());
        assert!(matches!(
            prepared.build(&numerators.with_zero_color_filter(true), &run).unwrap().into_amplitude(Default::default()),
            Err(GenerationError::Amplitude(error)) if matches!(*error, feynkit_amplitude::AmplitudeError::Empty)
        ));
    }

    #[test]
    fn zero_color_filter_progress_can_cancel_algebra() {
        let process = Process::new([] as [&str; 0], ["g"]);
        let run = Run::default().threads(1).progress(|update| {
            if update.stage == Stage::ZeroColor && update.completed == 1usize.into() {
                crate::GenerationControl::Cancel
            } else {
                crate::GenerationControl::Continue
            }
        });
        let result = process
            .amplitude(Loops::total(1..=1))
            .with_topology(
                TopologyRules::new()
                    .with_max_vertices(1)
                    .with_self_loops(true)
                    .with_zero_flow_edges(true),
            )
            .generate(Model::qcd().into(), &run)
            .unwrap()
            .build(
                &NumeratorOptions::default()
                    .with_naming(Naming::prefix("D"))
                    .with_zero_color_filter(true),
                &run,
            )
            .unwrap();
        assert!(!result.report().completed);
        assert_eq!(result.report().zero_numerator_count, 1usize.into());
        assert!(!result.diagrams().is_empty());
        result.validate_groups().unwrap();
    }

    #[test]
    fn single_thread_generation_stays_on_the_caller() {
        let caller = std::thread::current().id();
        let observed = Arc::new(Mutex::new(Vec::new()));
        let sink = observed.clone();
        let run = Run::default().threads(1).progress(move |progress| {
            if progress.completed > 0usize.into() {
                sink.lock()
                    .unwrap()
                    .push((progress.stage, std::thread::current().id()));
            }
            crate::GenerationControl::Continue
        });
        let generate = || {
            Process::new(["g"], ["g"])
                .amplitude(Loops::total(1..=1))
                .generate(Model::yang_mills().into(), &run)
                .unwrap()
                .build(
                    &NumeratorOptions::default().with_zero_color_filter(true),
                    &run,
                )
                .unwrap()
                .group(&Grouping::Identical(GraphGroupingOptions::default()), &run)
                .unwrap()
        };
        let first = generate();
        let second = generate();
        assert!(first.report().completed && second.report().completed);
        first.validate_groups().unwrap();
        second.validate_groups().unwrap();
        assert!(!first.is_empty());
        assert_eq!(first.len(), second.len());
        for (left, right) in first.iter().zip(second.iter()) {
            assert_eq!(left.to_dot().unwrap(), right.to_dot().unwrap());
        }
        let observed = observed.lock().unwrap();
        for stage in [
            Stage::Skeletons,
            Stage::Dressing,
            Stage::Prepare,
            Stage::Numerators,
            Stage::ZeroColor,
            Stage::GroupingPreparation,
        ] {
            assert!(
                observed.iter().any(|(actual, _)| *actual == stage),
                "missing {stage:?}"
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
            let run = Run::default().threads(threads).progress(move |progress| {
                if progress.stage == Stage::ZeroColor && progress.completed > 0usize.into() {
                    sink.lock().unwrap().push(rayon::current_num_threads());
                }
                crate::GenerationControl::Continue
            });
            let generated = Process::new(["g"], ["g"])
                .amplitude(Loops::total(1..=1))
                .with_topology(TopologyRules::new().with_self_loops(false))
                .generate(Model::yang_mills().into(), &run)
                .unwrap()
                .build(
                    &NumeratorOptions::default()
                        .with_naming(Naming::prefix("D"))
                        .with_zero_color_filter(true),
                    &run,
                )
                .unwrap();
            assert!(generated.report().completed);
            let observed = observed.lock().unwrap();
            assert!(!observed.is_empty());
            assert!(observed.iter().all(|&actual| actual == threads));
        }
    }

    /// Symbolica orders the factors of a product by symbol creation, so the
    /// printed numerators and overall factors must not depend on which diagram
    /// first meets a coupling or a factor head. Each run is a fresh process of
    /// this test binary: one loop at 1, 8 and 32 threads, and at 1 thread after
    /// the Born diagrams, which meet the external-ordering sign and the tree
    /// couplings before any loop diagram does.
    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn printed_numerators_and_overall_factors_do_not_depend_on_generation_order() {
        const RUN: &str = "FEYNKIT_GENERATOR_PRINT_TEST_RUN";
        if let Ok(run) = env::var(RUN) {
            let model = Arc::new(standard_model());
            let generate = |loops, threads| {
                let request = Request::amplitude(
                    Process::new(["e-", "e+"], ["mu-", "mu+"]),
                    crate::Loops::total(loops..=loops),
                )
                .with_topology(TopologyRules::default().with_self_loops(false));
                let run = Run::default().threads(threads);
                request
                    .generate(model.clone(), &run)
                    .unwrap()
                    .build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                    .unwrap()
            };
            let threads = match run.as_str() {
                "after-born" => {
                    generate(0, 1);
                    1
                }
                threads => threads.parse().unwrap(),
            };
            for diagram in generate(1, threads).iter() {
                println!(
                    "{RUN} {} {} {}",
                    diagram.name().expect("named build"),
                    NumeratorInterface::export(diagram.numerator().unwrap()).0,
                    diagram.overall_factor()
                );
            }
            return;
        }
        let child = |run: &str| {
            let output = Command::new(env::current_exe().unwrap())
                .args([
                    "--exact",
                    "generation::tests::printed_numerators_and_overall_factors_do_not_depend_on_generation_order",
                    "--nocapture",
                ])
                .env(RUN, run)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{run}: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            String::from_utf8(output.stdout)
                .unwrap()
                .lines()
                .filter_map(|line| line.split_once(RUN).map(|(_, record)| record.to_owned()))
                .collect::<Vec<_>>()
        };
        let serial = child("1");
        assert_eq!(serial.len(), 464);
        for run in ["8", "32", "after-born"] {
            let printed = child(run);
            assert_eq!(printed.len(), serial.len(), "{run}");
            for (serial, printed) in serial.iter().zip(&printed) {
                assert_eq!(printed, serial, "{run}");
            }
        }
    }

    #[test]
    fn numerator_progress_can_cancel_construction() {
        let generated = {
            let run = Run::default().threads(4).progress(|update| {
                if update.stage == Stage::Numerators && update.completed == 2usize.into() {
                    crate::GenerationControl::Cancel
                } else {
                    crate::GenerationControl::Continue
                }
            });
            let request = (Process::new(["g"], ["g"]))
                .clone()
                .amplitude(Loops::total(2..=2))
                .with_topology(TopologyRules::new().with_self_loops(false));
            request
                .generate((Model::yang_mills()).into(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        assert!(!generated.report().completed);
        assert!(!generated.diagrams().is_empty());
        assert!(
            generated.report().retained_count < generated.report().interaction_assignment_count
        );
        generated.validate_groups().unwrap();
    }

    #[test]
    fn cancellation_callback_returns_an_incomplete_result() {
        let generated = {
            let run = Run::default().cancellation_check(|| true);
            let request = (Process::new(["phi"], ["phi", "phi"]))
                .clone()
                .amplitude(Loops::total(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false));
            request
                .generate((scalar_model()).into(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();

        assert!(!generated.report().completed);
        assert!(generated.diagrams().is_empty());
    }

    #[test]
    fn vertex_signatures_present_ufo_particles_as_incoming_stubs() {
        let model = fermion_model();
        let options = photon_process()
            .amplitude(Loops::total(0..=0))
            .with_restrict(Restrict::new().only_vertices(["V"]));
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
            &(Process::new([1_i64, -1], Vec::<i64>::new()))
                .clone()
                .amplitude(Loops::total(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false)),
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
            &(Process::new([1_i64, -1], Vec::<i64>::new()))
                .clone()
                .amplitude(Loops::total(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false))
                .with_symmetrize(Symmetrize {
                    initial: true,
                    external_fermions: true,
                    ..Symmetrize::default()
                }),
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
            &(Process::new([1_i64, -1], Vec::<i64>::new()))
                .clone()
                .amplitude(Loops::total(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false)),
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
                &(Process::new([1000_i64], [1000_i64, 1000, 1000]))
                    .clone()
                    .cross_section(Loops::forward(3..=3))
                    .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                    .with_topology(TopologyRules::new().with_self_loops(false))
                    .with_cuts(
                        Cuts::new()
                            .with_blobs(1..=1)
                            .with_spectators(0..=0)
                            .with_symmetrize_left_right(enabled),
                    ),
            )
            .unwrap();
            let canonical = graphs
                .iter()
                .map(|graph| {
                    process
                        .canonicalize_numerator_graph(&graph.build(), &model, &Grouping::None)
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
            &(Process::new([1_i64, -1], [22_i64]))
                .clone()
                .cross_section(Loops::forward(0..=0))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false))
                .with_cuts(
                    Cuts::new()
                        .with_blobs(1..=1)
                        .with_spectators(0..=0)
                        .with_symmetrize_left_right(true),
                ),
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

        let grouping = Grouping::UpToScalar(GraphGroupingOptions::default());
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
            &(photon_process())
                .clone()
                .amplitude(Loops::total(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false)),
        )
        .unwrap();
        let grouping = Grouping::UpToSign(GraphGroupingOptions::default());

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

        let repaired = Generator::new(model.clone())
            .repair_cp_interaction_rules(graph.build())
            .unwrap();
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
        let repaired = Generator::new(model.clone())
            .repair_cp_interaction_rules(graph.build())
            .unwrap();
        assert_eq!(repaired.vertices().next().unwrap(), vertex(&model, "V_117"));

        let request =
            Process::new([-4_i64, 1, 251], Vec::<i64>::new()).amplitude(Loops::total(0..=0));
        let allow_original = request
            .clone()
            .with_restrict(Restrict::new().only_vertices(["V_82"]));
        let veto_repaired = request
            .clone()
            .with_restrict(Restrict::new().exclude_vertices(["V_117"]));
        let allow_repaired = request.with_restrict(Restrict::new().only_vertices(["V_117"]));
        let generator = Generator::new(model);

        assert!(
            !generator
                .passes_vertex_selection(&repaired, &allow_original)
                .unwrap(),
            "a CP repair must not bypass an allow-list"
        );
        assert!(
            !generator
                .passes_vertex_selection(&repaired, &veto_repaired)
                .unwrap(),
            "a CP repair must not bypass a veto-list"
        );
        assert!(
            generator
                .passes_vertex_selection(&repaired, &allow_repaired)
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
                &NumeratorOptions::default(),
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
    fn vertex_orders_take_each_order_maximum() {
        let model = fermion_model();
        let rule_orders = model.vertex_rule("V").unwrap().coupling_orders(&model);
        assert_eq!(
            rule_orders,
            BTreeMap::from([("QCD".to_owned(), 2), ("QED".to_owned(), 3)])
        );
    }

    fn admitted(expression: Atom) -> SymbolicTensor<PartialStructure> {
        NumeratorInterface::admit(expression, None).unwrap()
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
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(7), &legs[..3]);
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
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(7), &legs[..3]);
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
        assert!(
            generator.numerators.interactions.read().unwrap().is_empty(),
            "a rejected actual-owner binding is not cached"
        );
        assert!(
            generator
                .numerators
                .interaction_templates
                .read()
                .unwrap()
                .is_empty(),
            "an invalid intrinsic rule is not published"
        );
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
    fn ufo_lorentz_epsilon_matches_aloha_and_native_gamma5() {
        let mut legs = tensor_test_legs();
        for leg in &mut legs {
            leg.spin = 3;
        }
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(7), &legs);
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
        let legs: [NumeratorHalfEdge; 4] = std::array::from_fn(|i| NumeratorHalfEdge {
            edge: 10 * (i + 1),
            flow: Flow::Sink,
            hedge: 20 * (i + 1),
            spin: 3,
            color: 1,
        });
        let instantiation = NumeratorInstantiation::new(NumeratorOwner::Vertex(7), &legs);
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
        let options = photon_process()
            .amplitude(Loops::total(0..=0))
            .with_orders(Orders::new().exact("QCD", 6).exact("QED", 0));
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

        let generated = {
            let run = Run::default().threads(1);
            let request = (Process::new(["g"], ["g"]))
                .clone()
                .amplitude(Loops::total(1..=1))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(2),
                )
                .with_restrict(
                    Restrict::new().exclude_particles(
                        vec![-6, -4, -3, -2, -1, 1, 2, 3, 4, 6]
                            .into_iter()
                            .map(ParticleSelector::Pdg),
                    ),
                )
                .with_orders(Orders {
                    plain: BTreeMap::from([
                        ("QCD".to_owned(), (2, Some(2))),
                        ("QED".to_owned(), (0, Some(0))),
                    ])
                    .into_iter()
                    .map(|(name, (min, max))| (name, Range::Inclusive { min, max }))
                    .collect(),
                    ..Orders::new()
                });
            request.generate((model).into(), &run).and_then(|prepared| {
                prepared.build(
                    &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                    &run,
                )
            })
        }
        .unwrap();
        assert_eq!(generated.diagrams().len(), 3);
        assert_eq!(generated.report().retained_count, 3usize.into());
        assert_eq!(
            generated
                .diagrams()
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
        for diagram in generated.diagrams() {
            combined.push_str(
                &NumeratorInterface::export(diagram.numerator().unwrap())
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
        assert_eq!(
            diagram.numerator().unwrap().expression(),
            &(Atom::num(77) * test_atom("GC1"))
        );
        assert_eq!(
            diagram.loop_momentum_basis(),
            dot.diagram.loop_momentum_basis()
        );
        assert_eq!(diagram.projector().expression(), &Atom::one());
        diagram.validate().unwrap();
    }

    /// The external states of `diagram`'s projector as `(head, edge, half-edge)`.
    fn projector_states(diagram: &FeynmanDiagram) -> BTreeSet<(String, usize, usize)> {
        let mut states = BTreeSet::new();
        NumeratorInterface::export(diagram.projector())
            .0
            .as_view()
            .visitor(&mut |term| {
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
        assert_eq!(explicit.projector().expression(), &Atom::num(2));

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
        assert_eq!(diagram.numerator().unwrap().expression(), &Atom::one());
        let routed = diagram.loop_momentum_basis().clone();
        let id = diagram.id().unwrap();
        let diagram = diagram.apply_feynman_rules().unwrap();
        assert_eq!(
            diagram.numerator().unwrap().expression(),
            &test_atom("GC1^2")
        );
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
            numerator: admitted(Atom::num(5)),
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
        dummy.numerator = admitted(Atom::num(7));
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
                .numerator
                .expression(),
            &Atom::num(5)
        );
        assert_eq!(
            diagram
                .edges()
                .find(|(_, _, data)| data.is_dummy)
                .unwrap()
                .2
                .numerator
                .expression(),
            &Atom::num(7)
        );
    }

    #[test]
    fn apply_feynman_rules_matches_generated_local_factors() {
        let generator = Generator::new(standard_model());
        let generated = {
            let run = Run::default().threads(1);
            let request = (Process::new(["g"], ["g"]))
                .clone()
                .amplitude(Loops::total(1..=1))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(2),
                )
                .with_orders(Orders {
                    plain: BTreeMap::from([
                        ("QCD".to_owned(), (2, Some(2))),
                        ("QED".to_owned(), (0, Some(0))),
                    ])
                    .into_iter()
                    .map(|(name, (min, max))| (name, Range::Inclusive { min, max }))
                    .collect(),
                    ..Orders::new()
                });
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        assert!(!generated.diagrams().is_empty());
        for diagram in generated.diagrams() {
            let topology = diagram
                .clone()
                .with_numerator(admitted(Atom::one()))
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
                &NumeratorOptions::default(),
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
            &(Process::new([9000005_i64, 9000005], [9000005_i64, 9000005]))
                .clone()
                .amplitude(Loops::total(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false)),
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
        let generated = {
            let run = Run::default().threads(1);
            let request = (Process::new(["g"], ["g"]))
                .clone()
                .amplitude(Loops::total(1..=1))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_max_vertices(2),
                )
                .with_restrict(
                    Restrict::new().only_vertices(vec![VertexSelector::Name("V_35".to_owned())]),
                );
            request
                .generate(generator.model.clone(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        assert_eq!(generated.diagrams().len(), 1);
        let diagram = &generated.diagrams()[0];
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
                particle["spin"] = json!(2);
            }
        }
        for lorentz in definition["lorentz_structures"].as_array_mut().unwrap() {
            for spin in lorentz["spins"].as_array_mut().unwrap() {
                if spin.as_i64() == Some(-1) {
                    *spin = json!(2);
                }
            }
        }
        // This artificial spin/statistics model no longer describes the SM's
        // physical vector/Goldstone/ghost quartets used for covariant cuts.
        definition["covariant_cut_multiplets"] = json!({});
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
            let generated = {
                let run = Run::default().threads(1);
                let request = (Process::new(["g"], ["g"]))
                    .clone()
                    .amplitude(Loops::total(1..=1))
                    .with_topology(
                        TopologyRules::new()
                            .with_self_loops(false)
                            .with_max_vertices(2)
                            .with_fermion_loops(fermion_loops..=fermion_loops),
                    )
                    .with_restrict(
                        Restrict::new()
                            .only_vertices(vec![VertexSelector::Name("V_35".to_owned())]),
                    );
                request
                    .generate(generator.model.clone(), &run)
                    .and_then(|prepared| {
                        prepared.build(
                            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                            &run,
                        )
                    })
            }
            .unwrap();
            assert_eq!(generated.diagrams().len(), expected_diagrams);
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
                &NumeratorOptions::default(),
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
        let model = fermion_model();
        let process = photon_process();
        let amplitude = process.clone().amplitude(Loops::total(0..=0));
        let cross_section = process.clone().cross_section(Loops::forward(0..=0));
        let perturbative = Orders::new().perturbative("QED", 1);
        cross_section
            .clone()
            .with_orders(perturbative.clone())
            .validate_in(&model)
            .unwrap();
        assert!(
            matches!(amplitude.clone().with_orders(perturbative).validate_in(&model),
            Err(GenerationError::InvalidRequest(message)) if message.contains("require a cross section"))
        );
        let sewn = Cuts::new().with_sewn_tadpoles(true);
        cross_section
            .clone()
            .with_cuts(sewn.clone())
            .validate_in(&model)
            .unwrap();
        assert!(
            matches!(amplitude.clone().with_cuts(sewn).validate_in(&model),
            Err(GenerationError::InvalidRequest(message)) if message.contains("require a cross section"))
        );
        // A request has one explicit loop domain. Duplicate loop filters and cut-scoped
        // particle selectors are unrepresentable; reject a mismatched domain instead.
        assert!(
            matches!(Request::amplitude(process.clone(), Loops::forward(0..=1)).validate_in(&model),
            Err(GenerationError::InvalidRequest(message)) if message.contains("loop domain"))
        );
        assert!(
            matches!(Request::amplitude(process, Loops::total(std::ops::RangeInclusive::new(2, 1))).validate_in(&model),
            Err(GenerationError::InvalidRequest(message)) if message.contains("loop range is empty"))
        );
        for orders in [
            Orders::new().squared_exact("QED", 2),
            Orders::new().perturbative("QED", 1),
        ] {
            assert!(
                matches!(cross_section.clone().with_cuts(Cuts::new().with_cut_amplitude(CutAmplitude {
                orders, ..CutAmplitude::default()
            })).validate_in(&model), Err(GenerationError::InvalidRequest(message)) if message.contains("plain bounds only"))
            );
        }
    }

    #[test]
    fn grouping_options_are_validated() {
        let run = Run::default().threads(1);
        let built = || {
            Process::new(["phi", "phi"], ["phi", "phi"])
                .amplitude(Loops::total(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false))
                .generate(scalar_model().into(), &run)
                .unwrap()
                .build(
                    &NumeratorOptions::default().with_zero_color_filter(true),
                    &run,
                )
                .unwrap()
        };
        let invalid = Grouping::Identical(GraphGroupingOptions {
            number_of_numerical_samples: 0,
            ..GraphGroupingOptions::default()
        });
        assert!(matches!(
            built().group(&invalid, &run),
            Err(GenerationError::InvalidNumericalSampleCount)
        ));
        // Canonical comparison is an exact grouping mode without samples.
        let canonical_only = Grouping::Identical(GraphGroupingOptions {
            number_of_numerical_samples: 0,
            check_canonical_numerator: true,
            ..GraphGroupingOptions::default()
        });
        built().group(&canonical_only, &run).unwrap();
    }

    #[test]
    fn public_limits_report_arithmetic_overflow() {
        let model = fermion_model();
        let options = photon_process()
            .cross_section(Loops::forward(0..=0))
            .with_orders(
                Orders::new()
                    .perturbative("QED", usize::MAX)
                    .perturbative("QCD", 1),
            );
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
        let options = photon_process()
            .cross_section(Loops::forward(0..=0))
            .with_orders(Orders::new().perturbative("QED", 1));
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
            &(Process::new([-11_i64, 11], [5_i64, -5, 25]))
                .clone()
                .cross_section(Loops::forward(0..=0))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false)),
        )
        .unwrap();
        // These graphs contain the requested species, but neither has an
        // s-channel cut with the exact final-state multiset after allowing only
        // the configured unresolved QCD content.
        for (graph, qcd_order) in [(three_loop, 1), (four_loop, 2)] {
            let options = Process::new([-11_i64, 11], [5_i64, -5, 25])
                .cross_section(Loops::forward(0..=0))
                .with_orders(Orders::new().perturbative("QCD", qcd_order))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0));
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
            &(Process::new([22_i64, 22], [22_i64]))
                .clone()
                .cross_section(Loops::forward(0..=0))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false)),
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

        let options = Process::new([22_i64, 22], [22_i64])
            .cross_section(Loops::forward(0..=0))
            .with_cuts(
                Cuts::new()
                    .with_blobs(1..=1)
                    .with_spectators(0..=0)
                    .with_sewn_tadpoles(true),
            );
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
        assert!(
            process.passes_sewn_filter(
                &factorized.build(),
                &(Process::new([22_i64, 22], [22_i64]))
                    .clone()
                    .cross_section(Loops::forward(0..=0))
                    .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                    .with_topology(TopologyRules::new().with_self_loops(false))
                    .with_cuts(
                        Cuts::new()
                            .with_blobs(1..=1)
                            .with_spectators(0..=0)
                            .with_sewn_tadpoles(false)
                    )
            )
        );
    }

    /// The cross sections of `g g -> u u~` at one loop: their IDs, the sum of their evaluated
    /// overall factors, and the boxes among them as
    /// `(sewn copies on adjacent vertices, evaluated overall factor)`.
    fn gluon_fusion_boxes(symmetrize_initial: bool) -> (Vec<DiagramId>, Atom, Vec<(bool, Atom)>) {
        let process = Process::new(["g", "g"], ["u", "u~"]);
        let result = {
            let run = Run::default().threads(1);
            let request = (process)
                .clone()
                .cross_section(Loops::forward(1..=1))
                .with_restrict(Restrict::new().exclude_particles(["d", "c", "s", "t"]))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_topology(TopologyRules::new().with_self_loops(false))
                .with_symmetrize(Symmetrize {
                    initial: symmetrize_initial,
                    ..Symmetrize::default()
                });
            request
                .generate((Model::qcd()).into(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        let mut boxes = vec![];
        for diagram in result.diagrams() {
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
                .diagrams()
                .iter()
                .map(|diagram| diagram.id().unwrap())
                .collect(),
            Atom::add_many(result.diagrams().iter().map(evaluated_overall_factor)),
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
        request: &Request,
        symmetrize_initial: bool,
    ) -> Result<Vec<Atom>, GenerationError> {
        let mut factors = {
            let run = Run::default().threads(1);
            let request = request.clone().with_symmetrize(Symmetrize {
                initial: symmetrize_initial,
                ..Symmetrize::default()
            });
            request
                .generate((model.clone()).into(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }?
        .diagrams()
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
                .cross_section(Loops::forward(1..=1))
                .with_topology(TopologyRules::new().with_self_loops(false))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_restrict(Restrict::new().only_vertices(vertices.iter().copied()))
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
        let moller = Process::new(["e-", "e-"], ["e-", "e-"])
            .cross_section(Loops::forward(1..=1))
            .with_topology(TopologyRules::new().with_self_loops(false))
            .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0));
        for symmetrize_initial in [false, true] {
            assert_eq!(
                born_cross_section_factors(&Model::qed(), &moller, symmetrize_initial).unwrap(),
                factors(&[-1, 1])
            );
        }
        // a a -> t t~: |M_t|^2 and |M_u|^2 form one exchange orbit, M_t M_u* and M_u M_t* the
        // other.
        let photon_fusion = Process::new(["a", "a"], ["t", "t~"])
            .cross_section(Loops::forward(1..=1))
            .with_topology(TopologyRules::new().with_self_loops(false))
            .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
            .with_restrict(Restrict::new().only_particles(["a", "t"]));
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
                let all = {
                    let run = Run::default().threads(1);
                    let request = (process)
                        .clone()
                        .cross_section(Loops::forward(loops..=loops))
                        .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                        .with_topology(TopologyRules::new().with_self_loops(false))
                        .with_cuts(
                            Cuts::new()
                                .with_blobs(1..=1)
                                .with_spectators(0..=0)
                                .with_symmetrize_left_right(false),
                        )
                        .with_symmetrize(Symmetrize {
                            initial: symmetrize_initial,
                            ..Symmetrize::default()
                        });
                    request
                        .generate(generator.model.clone(), &run)
                        .and_then(|prepared| {
                            prepared.build(
                                &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                                &run,
                            )
                        })
                }
                .unwrap();
                assert!(all.report().completed);
                let mut expected = all
                    .diagrams()
                    .iter()
                    .filter(|diagram| selected(diagram))
                    .map(snapshot)
                    .collect::<Vec<_>>();
                assert!(!expected.is_empty(), "no QCD^2 diagrams at {loops} loops");
                assert!(
                    expected.len() < all.diagrams().len(),
                    "the budget must reject diagrams"
                );
                let filtered = {
                    let run = Run::default().threads(1);
                    let request = (process)
                        .clone()
                        .cross_section(Loops::forward(loops..=loops))
                        .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                        .with_topology(TopologyRules::new().with_self_loops(false))
                        .with_orders(Orders {
                            plain: BTreeMap::from([("QCD".to_owned(), (2, Some(2)))])
                                .into_iter()
                                .map(|(name, (min, max))| (name, Range::Inclusive { min, max }))
                                .collect(),
                            ..Orders::new()
                        })
                        .with_cuts(
                            Cuts::new()
                                .with_blobs(1..=1)
                                .with_spectators(0..=0)
                                .with_symmetrize_left_right(false),
                        )
                        .with_symmetrize(Symmetrize {
                            initial: symmetrize_initial,
                            ..Symmetrize::default()
                        });
                    request
                        .generate(generator.model.clone(), &run)
                        .and_then(|prepared| {
                            prepared.build(
                                &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                                &run,
                            )
                        })
                }
                .unwrap();
                assert!(filtered.report().completed);
                let mut actual = filtered.diagrams().iter().map(snapshot).collect::<Vec<_>>();
                expected.sort();
                actual.sort();
                assert_eq!(
                    actual, expected,
                    "{loops} loops, beam exchange {symmetrize_initial}"
                );

                // Compare each order class separately: signed totals across classes can cancel.
                let mut weights = BTreeMap::new();
                for diagram in all.diagrams() {
                    *weights.entry(orders(diagram)).or_insert_with(Atom::zero) +=
                        evaluated_overall_factor(diagram);
                }
                if symmetrize_initial {
                    assert_eq!(Some(weights), exact_weights, "{loops} loops");
                    // Two beams have only singleton and size-two orbits. Require the
                    // budgeted cases to exercise a size-two orbit, beyond fixed points.
                    nontrivial_orbit |= filtered.diagrams().iter().any(|diagram| {
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
        let generate = |process: Process, symmetrization: Symmetrize| {
            process
                .cross_section(Loops::forward(1..=1))
                .with_topology(TopologyRules::new().with_self_loops(false))
                .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                .with_symmetrize(symmetrization)
                .generate(model.clone().into(), &Run::default().threads(1))
        };
        assert!(matches!(
            generate(
                Process::new(["e-", "e+"], ["mu-", "mu+"]),
                Symmetrize::FINAL
            ),
            Err(GenerationError::CrossSectionFinalStateSymmetrization)
        ));
        // Three identical beams, or two among three incoming particles, have no exchange that
        // a Lorentz transformation realizes.
        for incoming in [["a", "a", "a"], ["a", "a", "e-"]] {
            assert!(matches!(
                generate(Process::new(incoming, ["e-", "e+"]), Symmetrize::INITIAL),
                Err(GenerationError::UnsupportedBeamSymmetrization { incoming: 3 })
            ));
        }
    }

    #[test]
    fn left_right_mirror_images_merge_into_one_diagram() {
        let process = Process::new(["e-", "e+"], ["mu-", "mu+"]);
        let generate = |left_right| {
            {
                let run = Run::default().threads(1);
                let request = (process)
                    .clone()
                    .cross_section(Loops::forward(2..=2))
                    .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
                    .with_topology(TopologyRules::new().with_self_loops(false))
                    .with_cuts(
                        Cuts::new()
                            .with_blobs(1..=1)
                            .with_spectators(0..=0)
                            .with_symmetrize_left_right(left_right),
                    );
                request
                    .generate((Model::qed()).into(), &run)
                    .and_then(|prepared| {
                        prepared.build(
                            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                            &run,
                        )
                    })
            }
            .unwrap()
            .diagrams()
            .to_vec()
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

        let result = {
            let run = Run::default().threads(1);
            let request = (Process::new(["phi"], Vec::<&str>::new()))
                .clone()
                .amplitude(Loops::total(2..=2))
                .with_topology(
                    TopologyRules::new()
                        .with_self_loops(false)
                        .with_bridges(Bridges::OnePIInLoops),
                );
            request
                .generate((scalar_model()).into(), &run)
                .and_then(|prepared| {
                    prepared.build(
                        &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                        &run,
                    )
                })
        }
        .unwrap();
        assert_eq!(result.diagrams().len(), 1);
        let diagram = &result.diagrams()[0];
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
        let classes = |request: Request| {
            ResolvedProcess::new(&model, request.process(), request.kind(), &request)
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
                process
                    .clone()
                    .amplitude(Loops::total(0..=0))
                    .with_symmetrize(Symmetrize::INITIAL)
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
                process
                    .amplitude(Loops::total(0..=0))
                    .with_symmetrize(Symmetrize::INITIAL | Symmetrize::EXTERNAL_FERMIONS)
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
                cross_section
                    .cross_section(Loops::forward(0..=0))
                    .with_symmetrize(Symmetrize::INITIAL)
                    .with_cuts(Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
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
                left_right.cross_section(Loops::forward(0..=0)).with_cuts(
                    Cuts::new()
                        .with_blobs(1..=1)
                        .with_spectators(0..=0)
                        .with_symmetrize_left_right(true)
                )
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
        fermion["propagator"] = json!("custom_f_prop");
        definition["propagators"]
            .as_array_mut()
            .unwrap()
            .push(json!({
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
                &NumeratorOptions::default(),
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

    #[test]
    fn routed_weights_match_all_symbolic_sign_components() {
        for mask in 0..32 {
            let signs = FermionSigns {
                internal_loop_negative: mask & 1 != 0,
                include_external_ordering: mask & 2 != 0,
                external_ordering_negative: mask & 4 != 0,
                include_antifermion_spin_sum: mask & 8 != 0,
                antifermion_spin_sum_negative: mask & 16 != 0,
            };
            let weight = DiagramWeight::new(signs, 19);
            let expected =
                feynkit_graph::expressions::evaluate_overall_factor(signs.factor(7, 19).as_view());
            // Weight evaluation does not need a graph. Read the independent exact components.
            assert_eq!(
                Atom::num(weight.sign) * Atom::num(weight.multiplicity) / Atom::num(7),
                expected,
                "mask {mask}"
            );
        }
    }

    #[test]
    fn routed_generation_matches_full_geometry_routing_and_grassmann_weights() {
        let topology = TopologyRules::new()
            .with_max_vertices(3)
            .with_self_loops(true);
        let cases = [
            (
                scalar_model(),
                Process::new(["phi"], ["phi"]).amplitude(Loops::total(0..=1)),
            ),
            (
                fermion_model(),
                Process::new([1_i64, 1], [1_i64, 1]).amplitude(Loops::total(0..=1)),
            ),
            (
                fermion_model(),
                Process::new([22_i64], [22_i64]).amplitude(Loops::total(0..=1)),
            ),
            (
                standard_model(),
                Process::new([-4_i64, 1, 251], std::iter::empty::<i64>())
                    .amplitude(Loops::total(0..=0)),
            ),
            (
                standard_model(),
                Process::new(["g"], ["g"])
                    .amplitude(Loops::total(1..=1))
                    .with_restrict(Restrict::new().only_vertices(["V_35"])),
            ),
        ];
        let run = Run::default().threads(1);
        let naming = Naming::prefix("D");
        let mut negative = 0;
        for (index, (model, request)) in cases.into_iter().enumerate() {
            let topology = if index == 3 {
                topology.clone().with_max_vertices(1)
            } else {
                topology.clone()
            };
            let prepared = request
                .with_topology(topology)
                .generate(Arc::new(model), &run)
                .unwrap();
            let full = prepared
                .build(
                    &NumeratorOptions::default().with_naming(naming.clone()),
                    &run,
                )
                .unwrap();
            let routed = prepared.routed(&naming).unwrap();
            assert!(!routed.is_empty());
            assert_eq!(full.len(), routed.len());
            for (full, routed) in full.iter().zip(&routed) {
                assert_eq!(full.id().unwrap(), routed.id().unwrap());
                assert_eq!(full.name(), routed.name());
                assert_eq!(full.loop_momentum_basis(), routed.loop_momentum_basis());
                assert_eq!(full.symmetry_factor(), routed.automorphism_order());
                assert_eq!(
                    numeric_factor(full),
                    Atom::num(routed.sign()) * Atom::num(routed.multiplicity())
                        / Atom::num(routed.automorphism_order())
                );
                assert_eq!(
                    full.vertices()
                        .map(|(id, data)| (id, data.name.as_str(), data.interaction))
                        .collect::<Vec<_>>(),
                    routed
                        .vertices()
                        .map(|(id, data)| (id, data.name, data.interaction))
                        .collect::<Vec<_>>()
                );
                assert_eq!(
                    full.edges()
                        .map(|(id, ends, data)| (
                            id,
                            ends,
                            data.particle,
                            data.directed,
                            data.external.as_ref(),
                            data.source_slot(),
                            data.target_slot()
                        ))
                        .collect::<Vec<_>>(),
                    routed
                        .edges()
                        .map(|(id, ends, data)| (
                            id,
                            ends,
                            data.particle,
                            data.directed,
                            data.external,
                            data.source_slot,
                            data.target_slot
                        ))
                        .collect::<Vec<_>>()
                );
                routed.validate().unwrap();
                negative += usize::from(routed.sign() < 0);
            }
        }
        assert!(
            negative > 0,
            "closed loops and exchanged external Grassmann lines carry actual minus signs"
        );
    }

    #[test]
    fn routed_generation_does_not_instantiate_local_rule_expressions() {
        let mut definition = serde_json::to_value(scalar_model()).unwrap();
        definition["lorentz_structures"][0]["structure"] = json!("P(1,99)");
        let model = Arc::new(Model::from_json(&definition.to_string()).unwrap());
        let run = Run::default().threads(1);
        let prepared = Process::new(["phi"], ["phi", "phi"])
            .amplitude(Loops::total(0..=0))
            .generate(model, &run)
            .unwrap();
        assert_eq!(prepared.routed(&Naming::None).unwrap().len(), 1);
        assert!(matches!(
            prepared.build(&NumeratorOptions::default(), &run),
            Err(GenerationError::InvalidNumeratorLeg { .. })
        ));
    }

    #[test]
    fn requested_projectors_name_process_external_legs() {
        let generator = Generator::new(scalar_model());
        let mut graph = Graph::new();
        let interaction = graph.add_node(ColoredNode::Interaction(vertex(&generator.model, "V1")));
        let external = graph.add_node(ColoredNode::External(ExternalNode {
            index: 0,
            state: ExternalState::Outgoing,
            particle: particle(&generator.model, 25),
        }));
        graph
            .add_edge(
                interaction,
                interaction,
                false,
                particle(&generator.model, 25),
            )
            .unwrap();
        graph
            .add_edge(interaction, external, false, particle(&generator.model, 25))
            .unwrap();
        // A requested projector names process external legs, not provisional edges: edge 0
        // of the generated graph is the self-loop here, external leg 0 is the carrier.
        let diagram = |projector: Atom| {
            generator.to_diagram(
                None,
                GenerationType::Amplitude,
                ColoredTopology {
                    graph: graph.build(),
                    symmetry: 1,
                    cut_partitions: vec![],
                    topology_threshold_candidates: vec![],
                },
                FermionSigns::default(),
                1,
                &NumeratorOptions::default().with_projector(projector),
            )
        };
        let projected = diagram(feynkit_graph::symbols::edge_index().call(0)).unwrap();
        let (carrier, _, _) = projected
            .edges()
            .find(|(_, _, data)| data.external.is_some())
            .unwrap();
        assert_eq!(
            carrier,
            EdgeId(0),
            "finalization places the external carrier first"
        );
        assert_eq!(
            NumeratorInterface::export(projected.projector()).0,
            feynkit_graph::symbols::edge_index().call(carrier.0)
        );
        assert!(matches!(
            diagram(feynkit_graph::symbols::edge_index().call(1)),
            Err(GenerationError::UnknownProjectorExternalIndex { index: 1 })
        ));
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
                .interaction_expression(rule, &instantiation, InteractionMomenta::Unreduced)
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
    fn interaction_momentum_reduction_public_api_retains_zero_interface() {
        let mut definition = serde_json::to_value(scalar_model()).unwrap();
        definition["lorentz_structures"][0]["structure"] =
            json!("P(idx(2,1),1)+P(idx(2,1),2)+P(idx(2,1),3)");
        let model = Model::from_json(&definition.to_string()).unwrap();
        let graph = FeynmanDiagram::from_dot(
            model,
            r#"digraph probe {
            ext [style=invis];
            ext -> a [particle="phi"];
            a -> ext [particle="phi"];
            a -> ext [particle="phi"];
        }"#,
        )
        .unwrap();
        let graph = graph.apply_feynman_rules().unwrap();
        assert_eq!(graph.vertices().count(), 1);
        for (_, vertex) in graph.vertices() {
            assert!(vertex.numerator.is_zero());
            assert_eq!(vertex.numerator.rank(), 1);
            assert_eq!(vertex.numerator.structure().logical_slots().len(), 1);
            assert!(
                vertex
                    .numerator
                    .explicit_indices()
                    .unwrap()
                    .next()
                    .is_some()
            );
        }
    }

    #[test]
    fn interaction_momentum_reduction_rejects_erased_invalid_indices() {
        for expression in [
            "P(idx(2,99),1)",
            "P(idx(2,99),1)+P(idx(2,99),2)+P(idx(2,99),3)",
        ] {
            let mut definition = serde_json::to_value(scalar_model()).unwrap();
            definition["lorentz_structures"][0]["structure"] = json!(expression);
            let model = Model::from_json(&definition.to_string()).unwrap();
            let rule = model.vertex_rule_id("V1").unwrap();
            let generator = Generator::new(model);
            let legs = tensor_test_legs();
            for owner in [4, 9] {
                let instantiation = NumeratorInstantiation {
                    owner: NumeratorOwner::Vertex(owner),
                    legs: &legs[..3],
                };
                for _ in 0..2 {
                    assert!(
                        matches!(
                            generator.interaction_numerator(rule, &instantiation),
                            Err(GenerationError::InvalidNumeratorLeg {
                                leg: 99,
                                legs: 3,
                                ..
                            })
                        ),
                        "{expression}"
                    );
                }
            }
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
    }

    #[test]
    fn interaction_momentum_reduction_rejects_heterogeneous_summands() {
        for expression in [
            "P(idx(2,1),1)+P(idx(2,1),2)+1",
            "P(idx(2,1),1)+P(idx(2,1),2)+P(idx(2,1),3)+1",
            "2*(P(idx(2,1),1)+P(idx(2,1),2)+P(idx(2,1),3))+1",
        ] {
            let mut definition = serde_json::to_value(scalar_model()).unwrap();
            definition["lorentz_structures"][0]["structure"] = json!(expression);
            let model = Model::from_json(&definition.to_string()).unwrap();
            let rule = model.vertex_rule_id("V1").unwrap();
            let generator = Generator::new(model);
            let legs = tensor_test_legs();
            for owner in [4, 9] {
                let instantiation = NumeratorInstantiation {
                    owner: NumeratorOwner::Vertex(owner),
                    legs: &legs[..3],
                };
                for _ in 0..2 {
                    assert!(
                        matches!(
                            generator.interaction_numerator(rule, &instantiation),
                            Err(GenerationError::Diagram(DiagramError::Tensor(
                                TensorInferenceError::Invalid(_)
                            )))
                        ),
                        "{expression}"
                    );
                }
            }
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
    }

    #[test]
    fn interaction_momentum_reduction_retains_complete_cancellation_metadata() {
        let _ = tensor_symbol!("UFO::generation_complete_color_left");
        let _ = tensor_symbol!("UFO::generation_complete_color_right");
        for color in [
            "1",
            // This closed pair is outside the Lorentz sum, so its construction
            // history must survive even though conservation erases all payload.
            "generation_complete_color_left(spenso::coad(8,dummy(1))) \
             *generation_complete_color_right(spenso::coad(8,dummy(1)))",
        ] {
            let mut definition = serde_json::to_value(scalar_model()).unwrap();
            definition["lorentz_structures"][0]["structure"] =
                json!("P(idx(2,1),1)+P(idx(2,1),2)+P(idx(2,1),3)");
            definition["vertex_rules"][0]["color_structures"][0] = json!(color);
            let model = Model::from_json(&definition.to_string()).unwrap();
            assert!(
                model.lorentz_structures()[0]
                    .reduced_by_momentum_conservation()
                    .is_zero()
            );
            let rule = model.vertex_rule_id("V1").unwrap();
            let generator = Generator::new(model);
            let mut legs = tensor_test_legs()[..3].to_vec();
            for leg in &mut legs {
                leg.color = 1;
            }
            let head = tensor_symbol!("generation_complete_momentum_sum_probe");
            let mut previous = None;
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
                assert!(zero.is_zero(), "{color}");
                assert_eq!(zero.rank(), 1, "{color}");
                assert!(zero.same_reduction_state(&cached));
                let public = AbstractIndex::Named(symbols::hedge_index().into(), legs[0].hedge, 2);
                assert_eq!(
                    zero.structure().logical_slots()[0].aind,
                    PartialIndex::Explicit(public)
                );
                let retained = cached.explicit_indices().unwrap().collect::<HashSet<_>>();
                assert!(retained.contains(&public));
                if color != "1" {
                    let dummy = AbstractIndex::Named(symbols::vertex_index().into(), owner, 1);
                    assert!(retained.contains(&dummy));
                    let probe = |index: AbstractIndex| {
                        SymbolicTensor::infer(
                            FunctionBuilder::new(head)
                                .add_arg(ColorAdjoint {}.new_rep(8).to_symbolic([index.to_atom()]))
                                .finish(),
                        )
                        .unwrap()
                    };
                    for numerator in [&zero, &cached] {
                        assert!(matches!(
                            numerator.multiply(&probe(dummy)),
                            Err(TensorCompositionError::InvalidExplicitMultiplicity {
                                index, occurrences: 3,
                            }) if index == dummy
                        ));
                        let other =
                            AbstractIndex::Named(symbols::vertex_index().into(), owner + 1, 1);
                        let product = numerator.multiply(&probe(other)).unwrap();
                        assert!(product.is_zero());
                        assert_eq!(product.rank(), 2);
                    }
                }
                let templates = generator.numerators.interaction_templates.read().unwrap();
                assert_eq!(templates.len(), 1);
                let template = Arc::clone(templates.values().next().unwrap());
                let NumeratorTemplate::Admitted { tensor, labels } = &*template else {
                    panic!("complete cancellation must retain shared raw-rule admission");
                };
                assert!(tensor.is_zero());
                assert_eq!(tensor.rank(), 1);
                assert!(!labels.is_empty());
                if let Some(previous) = previous {
                    assert!(Arc::ptr_eq(&previous, &template));
                }
                previous = Some(template);
                for leg in &mut legs {
                    leg.edge += 5;
                    leg.hedge += 7;
                }
            }
            assert_eq!(generator.numerators.interactions.read().unwrap().len(), 2);
        }
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
    fn interaction_templates_preserve_rank_when_incident_momenta_cancel() {
        let mut definition = serde_json::to_value(scalar_model()).unwrap();
        // Keep a shortest Lorentz rule: P1 + P2 would reduce to -P3 before
        // aliasing. Equal-flow aliases cancel P1 - P2 after admission instead.
        definition["lorentz_structures"][0]["structure"] = json!("P(idx(2,1),1)-P(idx(2,1),2)");
        let model = Model::from_json(&definition.to_string()).unwrap();
        let rule = model.vertex_rule_id("V1").unwrap();
        let generator = Generator::new(model);
        let mut legs = tensor_test_legs()[..3].to_vec();
        for leg in &mut legs {
            leg.color = 1;
        }
        legs[1].edge = legs[0].edge;
        legs[1].flow = legs[0].flow;
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
                    .interaction_expression(rule, &instantiation, InteractionMomenta::Unreduced)
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
    fn interaction_templates_use_canonical_named_interfaces_and_actual_owner_values() {
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
                        .interaction_expression(rule, &instantiation, InteractionMomenta::Unreduced)
                        .unwrap(),
                    None,
                )
                .unwrap()
                .canonicalize_named_interface();
                assert_eq!(
                    actual.expression().as_view().get_data(),
                    expected.expression().as_view().get_data(),
                    "{expression}: {permutation:?}"
                );
                assert_eq!(actual.presentation_atom(), *actual.expression());
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
                    .interaction_expression(rule, &instantiation, InteractionMomenta::Unreduced)
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
                .unwrap()
                .canonicalize_named_interface();
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
    fn propagator_templates_use_canonical_named_interfaces_and_actual_owner_values() {
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
                .unwrap()
                .canonicalize_named_interface();
                assert_eq!(
                    actual.expression().as_view().get_data(),
                    expected.expression().as_view().get_data(),
                    "{expression}: {hedges:?}"
                );
                assert_eq!(actual.presentation_atom(), *actual.expression());
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
            let direct =
                NumeratorInterface::translate(&source, |atom| instantiation.bind_template(atom))
                    .unwrap();
            let template = NumeratorTemplate::admitted(source).unwrap();
            let mapped = template.instantiate(&instantiation).unwrap();
            assert_eq!(
                mapped.explicit_indices().unwrap().collect::<Vec<_>>(),
                direct.explicit_indices().unwrap().collect::<Vec<_>>(),
                "template and direct translation must retain identical index history",
            );
            assert_eq!(
                mapped.structure().logical_slots(),
                vec![target[1], target[0]]
            );
            assert_eq!(mapped.is_zero(), was_zero);
            SymbolicTensor::validate_interface(mapped.expression(), mapped.structure()).unwrap();
        }
    }

    #[test]
    fn zero_template_and_direct_transport_retain_hidden_index_history() {
        let index = AbstractIndex::Named(symbols::hedge_index().into(), 0, 2);
        let port = Minkowski {}
            .new_rep(4)
            .cast::<LibraryRep>()
            .slot::<AbstractIndex, _>(index)
            .to_atom();
        let source = NumeratorInterface::admit(
            FunctionBuilder::new(tensor_symbol!("generation_template_hidden_pair"))
                .add_arg(&port)
                .add_arg(&port)
                .finish(),
            None,
        )
        .unwrap()
        .with_algebra_result(Atom::Zero)
        .unwrap();
        assert_eq!(source.rank(), 0);
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(7),
            legs: &legs[..2],
        };
        let direct =
            NumeratorInterface::translate(&source, |atom| instantiation.bind_template(atom))
                .unwrap();
        let mapped = NumeratorTemplate::admitted(source)
            .unwrap()
            .instantiate(&instantiation)
            .unwrap();
        let expected = vec![AbstractIndex::Named(
            symbols::hedge_index().into(),
            legs[0].hedge,
            2,
        )];
        assert_eq!(
            direct.explicit_indices().unwrap().collect::<Vec<_>>(),
            expected
        );
        assert_eq!(
            mapped.explicit_indices().unwrap().collect::<Vec<_>>(),
            expected
        );
        assert!(mapped.is_zero());
    }
}
