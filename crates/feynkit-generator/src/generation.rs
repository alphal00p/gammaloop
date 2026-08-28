use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    fmt,
    hash::{Hash, Hasher},
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    },
};

use feynkit_graph::{
    DiagramCut, DiagramCutSide, DiagramEdge, DiagramEndpoint, DiagramError, DiagramHalfEdge,
    DiagramId, DiagramVertex, EdgeId, ExternalState, FeynmanDiagram, VertexSlot,
};
use feynkit_model::{
    Model, ModelError, ParameterId, Particle, ParticleId, VertexRule, VertexRuleId,
};
use idenso::{
    color::CS,
    dirac::{AGS, PS},
    epsilon::EPSILON_SYMBOL,
    representations::{Bispinor, ColorAdjoint, ColorFundamental, ColorSextet},
};
use linnet::half_edge::involution::Flow;
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike};
use linnet::half_edge::{HedgeGraph, NodeIndex};
use rayon::{ThreadPoolBuildError, ThreadPoolBuilder, prelude::*};
use serde::{Deserialize, Serialize};
use spenso::{
    network::library::symbolic::ETS,
    structure::representation::{Minkowski, RepName},
};
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
    Process, ProcessError, SelectorError, SelfEnergyFilterOptions, SewnFilterOptions,
    SnailFilterOptions, TadpoleFilterOptions, VertexSelector, grouping, momentum_symbol,
};

#[derive(Debug, Error)]
pub enum GenerationError {
    #[error(transparent)]
    Process(#[from] ProcessError),
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

#[derive(Debug, Clone, Copy)]
struct NumeratorHalfEdge {
    edge: usize,
    flow: Flow,
    spin: i64,
    color: i64,
}

#[derive(Debug, Clone, Copy)]
enum NumeratorOwner {
    Vertex(usize),
    Edge(usize),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NumeratorSector {
    Spin,
    Color,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ColorRepresentation {
    Fundamental,
    AntiFundamental,
    Sextet,
    AntiSextet,
    Adjoint,
}

impl ColorRepresentation {
    fn dual(self) -> Self {
        match self {
            Self::Fundamental => Self::AntiFundamental,
            Self::AntiFundamental => Self::Fundamental,
            Self::Sextet => Self::AntiSextet,
            Self::AntiSextet => Self::Sextet,
            Self::Adjoint => Self::Adjoint,
        }
    }

    fn index(self, index: Atom) -> Atom {
        match self {
            Self::Fundamental => ColorFundamental {}.new_rep(3).to_symbolic([index]),
            Self::AntiFundamental => ColorFundamental {}.dual().new_rep(3).to_symbolic([index]),
            Self::Sextet => ColorSextet {}.new_rep(6).to_symbolic([index]),
            Self::AntiSextet => ColorSextet {}.dual().new_rep(6).to_symbolic([index]),
            Self::Adjoint => ColorAdjoint {}.new_rep(8).to_symbolic([index]),
        }
    }
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
                "C" => Err(self.unsupported_tensor(term)),
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

        // P(mu, leg) has already become FeynKit::Momentum(edge, mu).
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
        let _ = (CS.t, CS.f, ETS.metric);
        let mut representations = HashMap::<String, ColorRepresentation>::new();
        for leg in self.legs {
            if let Some(representation) = Self::vertex_color_representation(leg.color) {
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
                "f" => self.exact_arguments(term, &arguments, 3).and_then(|args| {
                    for argument in args {
                        self.assign_color_representation(
                            &mut representations,
                            argument.as_view(),
                            ColorRepresentation::Adjoint,
                        )?;
                    }
                    Ok(())
                }),
                "d" | "Epsilon" | "EpsilonBar" | "T6" | "K6" | "K6Bar" => {
                    Err(self.unsupported_tensor(term))
                }
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
                    ETS.metric(left.index(args[0].clone()), right.index(args[1].clone()))
                }),
                "T" => self.exact_arguments(term, &arguments, 3).map(|args| {
                    FunctionBuilder::new(CS.t)
                        .add_arg(ColorRepresentation::Adjoint.index(args[0].clone()))
                        .add_arg(ColorRepresentation::Fundamental.index(args[1].clone()))
                        .add_arg(ColorRepresentation::AntiFundamental.index(args[2].clone()))
                        .finish()
                }),
                "f" => self.exact_arguments(term, &arguments, 3).map(|args| {
                    args.iter()
                        .fold(FunctionBuilder::new(CS.f), |builder, argument| {
                            builder.add_arg(ColorRepresentation::Adjoint.index(argument.clone()))
                        })
                        .finish()
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
        let Some(source_representation) = Self::vertex_color_representation(color) else {
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
            source_representation.index(self.index(source, 1)?),
            source_representation.dual().index(self.index(sink, 1)?),
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
                NumeratorOwner::Vertex(_) => symbol!("FeynKit::VertexDummy"),
                NumeratorOwner::Edge(_) => symbol!("FeynKit::EdgeDummy"),
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
                    && function.get_nargs() == 2
        )
    }

    fn vertex_color_representation(color: i64) -> Option<ColorRepresentation> {
        match color {
            3 => Some(ColorRepresentation::Fundamental),
            -3 => Some(ColorRepresentation::AntiFundamental),
            6 => Some(ColorRepresentation::Sextet),
            -6 => Some(ColorRepresentation::AntiSextet),
            8 => Some(ColorRepresentation::Adjoint),
            _ => None,
        }
    }

    fn assign_color_representation(
        &self,
        representations: &mut HashMap<String, ColorRepresentation>,
        index: AtomView<'_>,
        representation: ColorRepresentation,
    ) -> Result<(), GenerationError> {
        let key = Self::index_key(index);
        match representations.insert(key, representation) {
            Some(previous) if previous != representation => Err(self.invalid_tensor(
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
        let options = options.resolve_selectors(&self.model)?;
        let options = &options;
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
            topology.graph = resolved.canonicalize_numerator_graph(
                &topology.graph,
                &self.model,
                &options.numerator_grouping,
            )?;
            // Fermion-flow normalization determines the physical sign, but
            // the canonical comparison graph remains the finalized graph.
            // Reorienting it here would change the momentum convention after
            // the reduced skeleton was selected and make mirror numerators
            // differ.  This matches the legacy ordering: canonicalize first,
            // instantiate on that representative, and use normalized flow as
            // sign metadata only.
            let normalized = resolved.normalize_fermion_flows(&topology.graph, &self.model)?;
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
                process.generation_type(),
                topology,
                fermion_loop_count,
                signs,
                options,
            )?);
        }
        if !options.forced_cuts.is_empty() {
            let mut matched = vec![false; options.forced_cuts.len()];
            let mut retained = Vec::new();
            for diagram in diagrams {
                let cuts = diagram
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
                    retained.push(diagram.with_cuts(cuts)?);
                }
            }
            diagrams = retained;
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
        let mut keyed_diagrams = diagrams
            .into_iter()
            .map(|diagram| Ok((diagram.canonical_key()?, diagram)))
            .collect::<Result<Vec<_>, GenerationError>>()?;
        keyed_diagrams.sort_by(|left, right| left.0.cmp(&right.0));
        let width = keyed_diagrams
            .len()
            .saturating_sub(1)
            .to_string()
            .len()
            .max(1);
        let mut diagrams = Vec::with_capacity(keyed_diagrams.len());
        for (index, (_, diagram)) in keyed_diagrams.into_iter().enumerate() {
            let name = format!("{}{index:0width$}", options.graph_prefix);
            let selected_by_id = options
                .selected_diagram_ids
                .as_ref()
                .is_some_and(|selected| selected.contains(&diagram.id()));
            let selected_by_name = options
                .selected_diagram_names
                .as_ref()
                .is_some_and(|selected| selected.contains(&name));
            let has_selection =
                options.selected_diagram_ids.is_some() || options.selected_diagram_names.is_some();
            if (has_selection && !selected_by_id && !selected_by_name)
                || options.vetoed_diagram_ids.contains(&diagram.id())
                || options.vetoed_diagram_names.contains(&name)
            {
                continue;
            }
            let diagram = diagram.with_name(&name);
            diagrams.push(
                if let Some(edges) = options.loop_momentum_bases.get(&diagram.id()) {
                    diagram.with_loop_momentum_edges(edges)?
                } else if let Some(edges) = options.named_loop_momentum_bases.get(&name) {
                    diagram.with_loop_momentum_edges(edges)?
                } else {
                    diagram
                },
            );
        }
        let grouped = if options.cancellation_requested() {
            completed = false;
            grouping::group_diagrams(
                diagrams,
                &self.model,
                &NumeratorGrouping::None,
                process.symmetrizes_left_right(),
            )?
        } else {
            grouping::group_diagrams(
                diagrams,
                &self.model,
                &options.numerator_grouping,
                process.symmetrizes_left_right(),
            )?
        };
        if options.cancellation_requested() {
            completed = false;
        }

        let result = GenerationResult {
            report: GenerationReport {
                topology_count,
                interaction_assignment_count,
                retained_count: grouped.diagrams.len(),
                zero_numerator_count: grouped.zero_numerator_count,
                completed,
            },
            diagrams: grouped.diagrams,
            groups: grouped.groups,
        };
        result.validate_groups()?;
        Ok(result)
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
            if !self.model.particle_by_id(edge.data.particle)?.is_fermion() {
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
        fermion_loop_count: usize,
        fermion_signs: FermionSigns,
        options: &GenerationOptions,
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
        let overall_factor = Atom::parse(
            &overall_factor,
            "feynkit_generator_factor",
            ParseSettings::default(),
        )
        .map_err(|message| GenerationError::NumeratorParse {
            owner: "diagram overall factor".to_owned(),
            expression: overall_factor,
            message,
        })?;
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
        let projector = match &options.projector {
            Some(projector) => projector.clone(),
            None => self.external_state_projector(&topology.graph)?,
        };
        let mut builder = FeynmanDiagram::builder(Arc::clone(&self.model), name)
            .symmetry_factor(topology.symmetry)
            .overall_factor(overall_factor)
            .numerator_prefactor(options.numerator_prefactor.clone())
            .projector(projector)
            .cuts(topology.cut_partitions.clone());
        let mut numerator_factors = Vec::new();
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
                    DiagramVertex::external_in_connection(
                        format!("ext{}", external.index),
                        external.index,
                        external.state,
                        connection,
                    )
                }
                ColoredNode::Interaction(rule_id) => {
                    let rule = self.model.vertex_rule_by_id(*rule_id)?;
                    let legs = self.interaction_half_edges(&topology.graph, index, rule)?;
                    let instantiation = NumeratorInstantiation {
                        owner: NumeratorOwner::Vertex(index),
                        legs: &legs,
                    };
                    for (slot, leg) in legs.iter().enumerate() {
                        endpoint_slots.insert((leg.edge, leg.flow), VertexSlot(slot));
                    }
                    let numerator = self.interaction_numerator(rule, &instantiation)?;
                    if numerator != Atom::one() {
                        numerator_factors.push(numerator.clone());
                    }
                    let mut vertex = DiagramVertex::interaction(format!("v{index}"), *rule_id);
                    vertex.numerator = numerator;
                    vertex
                }
            };
            builder.add_vertex(vertex);
        }
        for (edge_index, edge) in topology.graph.edges().iter().enumerate() {
            let particle = self.model.particle_by_id(edge.data.particle)?;
            let propagator_particle = if particle.is_antiparticle() {
                self.model.antiparticle(particle)?
            } else {
                particle
            };
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
                        spin: particle.spin,
                        color: particle.color,
                    },
                    NumeratorHalfEdge {
                        edge: edge_index,
                        flow: Flow::Sink,
                        spin: particle.spin,
                        color: particle.color,
                    },
                ];
                let instantiation = NumeratorInstantiation {
                    owner: NumeratorOwner::Edge(edge_index),
                    legs: &legs,
                };
                let propagator_numerator = self.propagator_numerator(propagator_particle)?;
                let spin =
                    instantiation.instantiate(&propagator_numerator, NumeratorSector::Spin)?;
                let color = instantiation.edge_color_identity(particle.color)?;
                match color {
                    None => spin,
                    Some(color) => spin * color,
                }
            } else {
                Atom::one()
            };
            if edge_numerator != Atom::one() {
                numerator_factors.push(edge_numerator.clone());
            }
            let mut diagram_edge = DiagramEdge::new(edge.data.particle, edge.directed);
            diagram_edge.numerator = edge_numerator;
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
        let numerator = numerator_factors
            .into_iter()
            .fold(Atom::one(), |product, factor| product * factor);
        let diagram = builder.numerator(numerator).build()?;
        diagram.validate()?;
        Ok(diagram)
    }

    fn external_state_projector(
        &self,
        graph: &Graph<ColoredNode, EdgeColor>,
    ) -> Result<Atom, GenerationError> {
        let mut projector = Atom::one();
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
            let index_symbol = match flow {
                Flow::Source => symbol!("FeynKit::SourceIndex"),
                Flow::Sink => symbol!("FeynKit::SinkIndex"),
            };
            let edge_index = i64::try_from(edge_id)
                .map_err(|_| GenerationError::ArithmeticOverflow("an external edge identifier"))?;
            let index = index_symbol.call((edge_index, 1_i64));
            let particle = self.model.particle_by_id(external.particle)?;
            let (head, index) = match particle.spin {
                2 => {
                    let head = match (external.state, particle.is_antiparticle()) {
                        (ExternalState::Incoming, false) => PS.u,
                        (ExternalState::Incoming, true) => PS.vbar,
                        (ExternalState::Outgoing, false) => PS.ubar,
                        (ExternalState::Outgoing, true) => PS.v,
                    };
                    (head, Bispinor {}.new_rep(4).to_symbolic([index]))
                }
                3 => {
                    let head = match external.state {
                        ExternalState::Incoming => PS.eps,
                        ExternalState::Outgoing => PS.ebar,
                    };
                    (head, Minkowski {}.new_rep(4).to_symbolic([index]))
                }
                _ => continue,
            };
            projector *= FunctionBuilder::new(head)
                .add_arg(edge_index)
                .add_arg(index)
                .finish();
        }
        Ok(projector)
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
            let particle = self.model.particle_by_id(edge.data.particle)?;
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
                        spin: 0,
                        color: 1,
                    },
                    base,
                    edge.directed.then_some(!particle.is_antiparticle()),
                ));
            }
            if edge.vertices.1 == vertex {
                available.push((
                    NumeratorHalfEdge {
                        edge: *edge_id,
                        flow: Flow::Sink,
                        spin: 0,
                        color: 1,
                    },
                    base,
                    edge.directed.then_some(particle.is_antiparticle()),
                ));
            }
        }

        let mut ordered = Vec::with_capacity(rule.particles.len());
        for (leg, particle_id) in rule.particles.iter().enumerate() {
            let particle = self.model.particle_by_id(*particle_id)?;
            let expected = ResolvedProcess::vertex_half_edge(&self.model, particle)?;
            let Some(position) = available.iter().position(|(_, particle, direction)| {
                *particle == expected.data.particle && *direction == expected.direction
            }) else {
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

    fn propagator_numerator(&self, particle: &Particle) -> Result<Atom, GenerationError> {
        if let Some(propagator) = particle.propagator {
            return Ok(self.model.propagator_by_id(propagator)?.numerator.clone());
        }
        let antiparticle = self.model.antiparticle(particle)?;
        if let Some(propagator) = antiparticle.propagator {
            return Ok(self.model.propagator_by_id(propagator)?.numerator.clone());
        }
        let particle_id = self.model.particle_id(&particle.name)?;
        Ok(self
            .model
            .propagators()
            .iter()
            .find(|propagator| {
                propagator.particle == particle_id || propagator.particle == particle.antiparticle
            })
            .map_or_else(Atom::one, |propagator| propagator.numerator.clone()))
    }

    fn interaction_numerator(
        &self,
        rule: &VertexRule,
        instantiation: &NumeratorInstantiation<'_>,
    ) -> Result<Atom, GenerationError> {
        let mut terms = Vec::new();
        for (color_index, color) in rule.color_structures.iter().enumerate() {
            for (lorentz_index, lorentz_name) in rule.lorentz_structures.iter().enumerate() {
                if let Some(coupling) = rule
                    .couplings
                    .get(color_index)
                    .and_then(|row| row.get(lorentz_index))
                    .and_then(Option::as_ref)
                {
                    let coupling = &self.model.coupling_by_id(*coupling)?.expression;
                    let lorentz = &self.model.lorentz_structure_by_id(*lorentz_name)?.structure;
                    let color = instantiation.instantiate(color, NumeratorSector::Color)?;
                    let lorentz = instantiation.instantiate(lorentz, NumeratorSector::Spin)?;
                    terms.push(coupling.clone() * color * lorentz);
                }
            }
        }
        if terms.is_empty() {
            Ok(Atom::one())
        } else {
            Ok(terms.into_iter().fold(Atom::new(), |sum, term| sum + term))
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct EdgeColor {
    particle: ParticleId,
    direction: Option<bool>,
}

impl fmt::Display for EdgeColor {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "particle#{}", self.particle.index())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct ExternalNode {
    index: usize,
    state: ExternalState,
    particle: ParticleId,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum ExternalCanonicalClass {
    Exact { state: ExternalState, index: usize },
    WithinState(ExternalState),
}

#[derive(Debug, Clone, Default)]
struct NodeColor {
    external: Option<ExternalNode>,
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

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum NumeratorCanonicalNode {
    Internal,
    External(ExternalNode),
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum NumeratorCanonicalEdge {
    InternalMassAndSpin {
        mass: ParameterId,
        spin: i64,
    },
    InternalParticle(ParticleId),
    External {
        particle: ParticleId,
        directed: bool,
    },
}

type InteractionAssignment = (Graph<ColoredNode, EdgeColor>, usize);

impl fmt::Display for ColoredNode {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::External(external) => write!(formatter, "ext{}", external.index),
            Self::Interaction(rule) => write!(formatter, "vertex#{}", rule.index()),
        }
    }
}

struct ColoredTopology {
    graph: Graph<ColoredNode, EdgeColor>,
    symmetry: u64,
    multiplicity: usize,
    cut_partitions: Vec<DiagramCut>,
}

#[derive(Debug, Clone, Copy, Default)]
struct FermionSigns {
    include_external_ordering: bool,
    external_ordering_negative: bool,
    include_antifermion_spin_sum: bool,
    antifermion_spin_sum_negative: bool,
}

struct NormalizedFermionFlow {
    #[cfg_attr(not(test), allow(dead_code))]
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
    fn new(model: &Model, process: &Process) -> Result<Self, GenerationError> {
        let resolve = |selector: &ParticleSelector| -> Result<Particle, GenerationError> {
            Ok(model.particle_by_id(selector.resolve(model)?)?.clone())
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
        let symmetrized_within_state = match state {
            ExternalState::Incoming => self.symmetrize_initial,
            ExternalState::Outgoing => self.symmetrize_final,
        } && (self.generation_type == GenerationType::CrossSection
            || !particle.is_fermion()
            || self.symmetrize_external_fermions);
        // Left/right symmetry is one global transformation of a sewn graph,
        // not an independent automorphism of every external pair. It is
        // applied after generation when canonical diagram keys are built.
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
            NumeratorGrouping::None | NumeratorGrouping::OnlyDetectZeroes => None,
        };
        let include_left_right =
            self.generation_type == GenerationType::CrossSection && self.symmetrize_left_right;
        let exact_internal_particles =
            grouping_options.is_some_and(|options| !options.differentiate_particle_masses_only);

        let mut candidates = Vec::new();
        let mut remappings = vec![BTreeMap::new()];
        if include_left_right {
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
                    "left/right canonicalization expected {} sewn external nodes, found {}",
                    2 * incoming_count,
                    external_nodes.len()
                )));
            }
            remappings.push(
                (0..incoming_count)
                    .flat_map(|index| {
                        [
                            (
                                external_nodes[&index],
                                external_nodes[&(incoming_count + index)],
                            ),
                            (
                                external_nodes[&(incoming_count + index)],
                                external_nodes[&index],
                            ),
                        ]
                    })
                    .collect(),
            );
        }
        for (swapped, remapping) in remappings.into_iter().enumerate() {
            let mut skeleton = Graph::new();
            for node in graph.nodes() {
                skeleton.add_node(match &node.data {
                    ColoredNode::External(external) => {
                        NumeratorCanonicalNode::External(external.clone())
                    }
                    ColoredNode::Interaction(_) => NumeratorCanonicalNode::Internal,
                });
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
                // into the interaction graph.  Internal direction is omitted
                // from the reduced skeleton but restored below.
                if matches!(
                    graph.nodes()[edge.vertices.1].data,
                    ColoredNode::External(_)
                ) {
                    vertices = (vertices.1, vertices.0);
                }
                let color = if external {
                    NumeratorCanonicalEdge::External {
                        particle: edge.data.particle,
                        directed: edge.directed,
                    }
                } else if exact_internal_particles {
                    NumeratorCanonicalEdge::InternalParticle(edge.data.particle)
                } else {
                    let particle = model.particle_by_id(edge.data.particle)?;
                    NumeratorCanonicalEdge::InternalMassAndSpin {
                        mass: particle.mass,
                        spin: particle.spin,
                    }
                };
                skeleton
                    .add_edge(vertices.0, vertices.1, external && edge.directed, color)
                    .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
            }
            candidates.push((swapped == 1, remapping, skeleton.canonize()));
        }
        candidates.sort_by(|left, right| left.2.graph.cmp(&right.2.graph));
        let mut representatives = Vec::new();
        for (left_right_swapped, remapping, canonical) in candidates {
            let mut canonical_to_input = vec![0; graph.nodes().len()];
            for (input, canonical_position) in canonical.vertex_map.iter().copied().enumerate() {
                canonical_to_input[canonical_position] =
                    remapping.get(&input).copied().unwrap_or(input);
            }
            let mut input_to_canonical = vec![0; graph.nodes().len()];
            for (canonical_position, input) in canonical_to_input.iter().copied().enumerate() {
                input_to_canonical[input] = canonical_position;
            }

            let mut result = Graph::new();
            for input in canonical_to_input.iter().copied() {
                // The selected global remapping is an involution. Applying it
                // a second time keeps each physical external label attached to
                // its original process state while adopting the selected
                // topology.
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
                    if edge.directed && left_right_swapped {
                        data.particle = model.particle_by_id(data.particle)?.antiparticle;
                    }
                    Ok((vertices, edge.directed, data))
                })
                .collect::<Result<Vec<_>, GenerationError>>()?;
            edges.sort_by_key(|(vertices, directed, data)| {
                (vertices.0, vertices.1, data.particle, *directed)
            });
            for (vertices, directed, data) in edges {
                result
                    .add_edge(vertices.0, vertices.1, directed, data)
                    .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
            }
            if include_left_right {
                // After the external involution is applied twice, exact
                // process labels are restored.  In the base/anti ParticleId
                // representation this leaves one residual CP ambiguity on
                // the internal charged flows.  Add exactly its one global
                // conjugate; taking the minimum below must never conjugate a
                // subset of those edges independently.
                let mut conjugated = Graph::new();
                for node in result.nodes() {
                    conjugated.add_node(node.data.clone());
                }
                for edge in result.edges() {
                    let external = matches!(
                        result.nodes()[edge.vertices.0].data,
                        ColoredNode::External(_)
                    ) || matches!(
                        result.nodes()[edge.vertices.1].data,
                        ColoredNode::External(_)
                    );
                    let mut data = edge.data;
                    if edge.directed && !external {
                        data.particle = model.particle_by_id(data.particle)?.antiparticle;
                    }
                    conjugated
                        .add_edge(edge.vertices.0, edge.vertices.1, edge.directed, data)
                        .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
                }
                representatives.push(conjugated);
            }
            representatives.push(result);
        }
        representatives
            .into_iter()
            .min()
            .ok_or_else(|| GenerationError::GraphConstruction("no canonical graph".to_owned()))
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
        // Pairing of the external fermion flows. The keys are sorted pairs of
        // all-incoming PDGs, and values are connected external leg IDs in the
        // order of that key.
        let mut pairings: BTreeMap<(i64, i64), Vec<(usize, usize)>> = BTreeMap::new();
        let mut line_pairings = Vec::new();

        // First fix flows connected to external legs, then all internal flows.
        for (node_id, external) in external_nodes.values() {
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
            // Ghost chains are oriented like fermion chains but do not enter
            // the external-fermion permutation sign. Physical processes do
            // not normally expose external ghosts, yet keeping this split
            // makes the flow normalization faithful for complete UFO models.
            if !model.particle_by_id(external.particle)?.is_fermion() {
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
            let particle = model.particle_by_id(graph.edges()[edge_id].data.particle)?;
            if (particle.is_fermion() || particle.is_ghost()) && !visited_edges[edge_id] {
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
            let mut data = edge.data;
            if edge.directed && (left, right) != edge.vertices {
                data.particle = model.particle_by_id(data.particle)?.antiparticle;
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
    ) -> Result<Vec<DiagramCut>, GenerationError> {
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
            partitions.push(DiagramCut {
                cut: stable_half_edges(&cut.left),
                left: summarize_side(&left)?,
                right: summarize_side(&right)?,
            });
        }
        Ok(partitions)
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
    use crate::{CancellationToken, GraphGroupingOptions, NumeratorGrouping};
    use feynkit_model::Model;

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
                spin: 2,
                color: -3,
            },
            NumeratorHalfEdge {
                edge: 20,
                flow: Flow::Sink,
                spin: 2,
                color: 3,
            },
            NumeratorHalfEdge {
                edge: 30,
                flow: Flow::Source,
                spin: 3,
                color: 8,
            },
            NumeratorHalfEdge {
                edge: 40,
                flow: Flow::Sink,
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

    fn standard_model() -> Model {
        Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap()
    }

    fn photon_process() -> Process {
        Process::amplitude([22_i64], [22_i64])
    }

    fn scalar_cross_section_process() -> Process {
        Process::cross_section(["phi", "phi"], ["phi", "phi"])
            .with_loop_count(1, 1)
            .unwrap()
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
                .all(|diagram| !diagram.numerator().is_zero())
        );
    }

    #[test]
    fn finalizes_external_fermion_projectors_unless_explicitly_overridden() {
        let generator = Generator::new(fermion_model());
        let process = Process::amplitude([1_i64, -1], [1_i64, -1])
            .with_loop_count(0, 0)
            .unwrap();
        let options = GenerationOptions::default().threads(1).max_vertices(2);
        let generated = generator.generate(&process, &options).unwrap();
        assert!(!generated.diagrams.is_empty());
        for diagram in &generated.diagrams {
            for wavefunction in [PS.u, PS.ubar, PS.v, PS.vbar] {
                assert_eq!(function_count(diagram.projector(), wavefunction), 1);
            }
            assert_eq!(
                function_count(diagram.projector(), symbol!("FeynKit::SourceIndex"))
                    + function_count(diagram.projector(), symbol!("FeynKit::SinkIndex")),
                4
            );
        }

        let explicit = generator
            .generate(&process, &options.projector(Atom::one()))
            .unwrap();
        assert!(
            explicit
                .diagrams
                .iter()
                .all(|diagram| diagram.projector() == &Atom::one())
        );
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
        let vector_projector = generator.external_state_projector(&vector_graph).unwrap();
        assert_eq!(function_count(&vector_projector, PS.eps), 1);
        assert_eq!(function_count(&vector_projector, PS.ebar), 1);
        assert_eq!(
            function_count(&vector_projector, symbol!("FeynKit::SinkIndex")),
            1
        );
        assert_eq!(
            function_count(&vector_projector, symbol!("FeynKit::SourceIndex")),
            1
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
            generator.external_state_projector(&scalar_graph).unwrap(),
            Atom::one()
        );
    }

    #[test]
    fn forced_cuts_retain_only_the_requested_typed_edge_set() {
        let generator = Generator::new(scalar_model());
        let process = scalar_cross_section_process();
        let options = GenerationOptions::default().threads(1).max_vertices(4);
        let generated = generator.generate(&process, &options).unwrap();
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
            .generate(&process, &options.clone().forced_cuts([requested.clone()]))
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
                .threads(1)
                .max_vertices(4)
                .forced_cuts([unknown]),
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
        let options = GenerationOptions::default().threads(1).max_vertices(4);
        let generated = generator.generate(&process, &options).unwrap();
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
        let process = Process::amplitude([foreign_particle.clone()], ["phi"]);
        assert!(matches!(
            Generator::new(target.clone()).generate(&process, &GenerationOptions::default()),
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
                Generator::new(target.clone())
                    .generate(&Process::amplitude(["phi"], ["phi", "phi"]), &options,),
                Err(GenerationError::Selector(
                    SelectorError::ModelMismatch { .. }
                ))
            ));
        }
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
        let process =
            ResolvedProcess::new(&model, &Process::amplitude([1_i64, -1], Vec::<i64>::new()))
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
        graph.add_edge(right, left, true, fermion).unwrap();
        graph.add_edge(right, antiparticle, true, fermion).unwrap();

        let normalized = process.normalize_fermion_flows(&graph, &model).unwrap();
        assert_eq!(
            normalized.graph.edges()[0].vertices,
            (particle_vertex, left)
        );
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
    fn ghost_chains_are_oriented_without_entering_external_fermion_signs() {
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
        let process =
            ResolvedProcess::new(&model, &Process::amplitude([1_i64, -1], Vec::<i64>::new()))
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
        assert!(!normalized.signs.external_ordering_negative);
    }

    #[test]
    fn global_left_right_canonicalization_preserves_charged_slash_momenta() {
        let model = fermion_model();
        let process = ResolvedProcess::new(
            &model,
            &Process::cross_section([1_i64, -1], [22_i64]).symmetrize_left_right(true),
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
        let particle_numerator = generator
            .propagator_numerator(model.particle_by_id(particle(&model, 1)).unwrap())
            .unwrap();
        let antiparticle_numerator = generator
            .propagator_numerator(model.particle_by_id(particle(&model, -1)).unwrap())
            .unwrap();
        assert_eq!(particle_numerator, antiparticle_numerator);
        let propagator = generator
            .propagator_numerator(
                model
                    .particle_by_id(canonical.edges()[internal_edge].data.particle)
                    .unwrap(),
            )
            .unwrap();
        let legs = [
            NumeratorHalfEdge {
                edge: internal_edge,
                flow: Flow::Source,
                spin: 2,
                color: 1,
            },
            NumeratorHalfEdge {
                edge: internal_edge,
                flow: Flow::Sink,
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
        assert_eq!(function_count(&localized, symbol!("FeynKit::Momentum")), 1);
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
    fn interaction_numerator_uses_model_expressions() {
        let model = fermion_model();
        let generator = Generator::new(model.clone());
        let legs = tensor_test_legs();
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(7),
            legs: &legs[..3],
        };
        let numerator = generator
            .interaction_numerator(model.vertex_rule("V").unwrap(), &instantiation)
            .unwrap();

        let numerator = numerator.to_plain_string();
        assert!(numerator.contains("spenso::gamma"));
        assert!(numerator.contains("spenso::t"));
        assert!(!numerator.contains("Gamma("));
        assert!(!numerator.contains("UFO::T("));
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
                "spenso::gamma(spenso::bis(4,FeynKit::SinkIndex(20,1)),spenso::bis(4,FeynKit::SourceIndex(10,1)),spenso::mink(4,FeynKit::SourceIndex(30,1)))"
            )
        );

        let metric = instantiation
            .instantiate(&test_atom("Metric(3,4)"), NumeratorSector::Spin)
            .unwrap();
        assert_eq!(
            metric,
            test_atom(
                "spenso::g(spenso::mink(4,FeynKit::SourceIndex(30,1)),spenso::mink(4,FeynKit::SinkIndex(40,1)))"
            )
        );

        let identity = instantiation
            .instantiate(&test_atom("Identity(1,2)"), NumeratorSector::Spin)
            .unwrap();
        assert_eq!(
            identity,
            test_atom(
                "spenso::g(spenso::bis(4,FeynKit::SourceIndex(10,1)),spenso::bis(4,FeynKit::SinkIndex(20,1)))"
            )
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
            "spenso::g(spenso::dind(spenso::cof(3,FeynKit::SourceIndex(10,1))),spenso::cof(3,FeynKit::SinkIndex(20,1)))*spenso::t(spenso::coad(8,FeynKit::SourceIndex(30,1)),spenso::cof(3,FeynKit::SinkIndex(20,1)),spenso::dind(spenso::cof(3,FeynKit::SourceIndex(10,1))))*spenso::f(spenso::coad(8,FeynKit::SourceIndex(30,1)),spenso::coad(8,FeynKit::VertexDummy(7,1)),spenso::coad(8,FeynKit::VertexDummy(7,2)))",
        );
        assert_eq!(tensors, expected);
    }

    #[test]
    fn expands_pslash_with_one_shared_collision_free_lorentz_index() {
        let legs = [
            NumeratorHalfEdge {
                edge: 5,
                flow: Flow::Source,
                spin: 2,
                color: 3,
            },
            NumeratorHalfEdge {
                edge: 5,
                flow: Flow::Sink,
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
                "spenso::gamma(spenso::bis(4,FeynKit::SinkIndex(5,1)),spenso::bis(4,FeynKit::SourceIndex(5,1)),spenso::mink(4,FeynKit::EdgeDummy(5,5)))*FeynKit::Momentum(5,spenso::mink(4,FeynKit::EdgeDummy(5,5)))*spenso::g(spenso::mink(4,FeynKit::EdgeDummy(5,4)),spenso::mink(4,FeynKit::EdgeDummy(5,4)))"
            )
        );
    }

    #[test]
    fn keeps_color_and_spin_identity_sectors_distinct() {
        let model = standard_model();
        let generator = Generator::new(model.clone());
        let legs = [
            NumeratorHalfEdge {
                edge: 1,
                flow: Flow::Source,
                spin: 2,
                color: -3,
            },
            NumeratorHalfEdge {
                edge: 2,
                flow: Flow::Sink,
                spin: 2,
                color: 3,
            },
            NumeratorHalfEdge {
                edge: 3,
                flow: Flow::Source,
                spin: 1,
                color: 1,
            },
        ];
        let instantiation = NumeratorInstantiation {
            owner: NumeratorOwner::Vertex(9),
            legs: &legs,
        };
        let numerator = generator
            .interaction_numerator(model.vertex_rule("V_78").unwrap(), &instantiation)
            .unwrap();

        let numerator = numerator.to_plain_string();
        assert!(numerator.contains("spenso::bis(4"));
        assert!(numerator.contains("spenso::cof(3"));
        assert!(numerator.contains("spenso::dind"));
        assert_eq!(numerator.matches("spenso::g(").count(), 2);
        assert!(!numerator.contains("UFO::Identity("));
    }

    #[test]
    fn generated_qcd_loop_numerators_are_spenso_ready() {
        let model = standard_model();
        let options = GenerationOptions::default()
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
                &Process::amplitude(["g"], ["g"])
                    .with_loop_count(1, 1)
                    .unwrap(),
                &options,
            )
            .unwrap();
        assert_eq!(generated.diagrams.len(), 3);

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

        let loop_count = generator.closed_fermion_loop_count(&graph).unwrap();
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
                },
                loop_count,
                FermionSigns::default(),
                &GenerationOptions::default(),
            )
            .unwrap();
        assert_eq!(
            diagram.overall_factor(),
            &factor_atom("InternalFermionLoopSign(-1)*CouplingsMultiplicity(2)/AutG(3)")
        );
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
                },
                0,
                FermionSigns {
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

        let invalid_scope = GenerationOptions::default().with_cut_amplitude_filter(
            GenerationFilter::ParticleVeto(vec![ParticleSelector::Pdg(1)]),
        );
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
        let classes = |process: Process| {
            ResolvedProcess::new(&model, &process)
                .unwrap()
                .external_edges(&model)
                .unwrap()
                .into_iter()
                .map(|(node, _)| node.external_class.unwrap())
                .collect::<Vec<_>>()
        };
        let process =
            Process::amplitude([1_i64, -1, 22], Vec::<i64>::new()).symmetrize_initial(true);
        assert_eq!(
            classes(process.clone()),
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
            classes(process.symmetrize_external_fermions(true)),
            vec![
                ExternalCanonicalClass::WithinState(ExternalState::Incoming),
                ExternalCanonicalClass::WithinState(ExternalState::Incoming),
                ExternalCanonicalClass::WithinState(ExternalState::Incoming),
            ]
        );

        let cross_section = Process::cross_section([1_i64, -1], [22_i64]).symmetrize_initial(true);
        assert_eq!(
            classes(cross_section),
            vec![
                ExternalCanonicalClass::WithinState(ExternalState::Incoming),
                ExternalCanonicalClass::WithinState(ExternalState::Incoming),
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

        let left_right = Process::cross_section([1_i64, -1], [22_i64]).symmetrize_left_right(true);
        assert_eq!(
            classes(left_right),
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
                },
                1,
                FermionSigns::default(),
                &GenerationOptions::default(),
            )
            .unwrap();
        let internal_numerator = diagram
            .edges()
            .find(|(id, _, _)| id.0 == 0)
            .unwrap()
            .2
            .numerator
            .to_plain_string();
        assert!(internal_numerator.contains("FeynKit::Momentum"));
        assert!(internal_numerator.contains("FeynKit::SourceIndex"));
        assert!(internal_numerator.contains("FeynKit::SinkIndex"));
        assert!(internal_numerator.contains("FeynKit::EdgeDummy"));
        assert!(!internal_numerator.contains("idx("));
        assert!(!internal_numerator.contains("Slash(P)"));

        let momentum = symbol!("FeynKit::Momentum");
        assert!(momentum.has_tag("spenso::tensor"));
        assert!(momentum.has_tag("spenso::rank1"));

        let external_numerator = &diagram
            .edges()
            .find(|(id, _, _)| id.0 == 1)
            .unwrap()
            .2
            .numerator;
        assert_eq!(external_numerator, &Atom::one());

        let vertex_numerator = diagram
            .vertex(feynkit_graph::VertexId(interaction))
            .unwrap()
            .numerator
            .to_plain_string();
        assert!(vertex_numerator.contains("FeynKit::SourceIndex"));
        assert!(vertex_numerator.contains("FeynKit::SinkIndex"));
    }
}
