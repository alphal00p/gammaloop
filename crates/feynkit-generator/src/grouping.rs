use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
    sync::atomic::{AtomicUsize, Ordering},
};

use feynkit_graph::{
    DiagramEndpoint, DiagramError, DiagramHalfEdge, EdgeId, ExternalLeg, ExternalState,
    FeynmanDiagram, VertexId,
};
use feynkit_model::{Model, ModelError, ParameterId, ParticleId};
use idenso::{
    IndexTooling,
    color::{CS, ColorSimplifier, ColorSimplifySettings},
    dirac::{AGS, PS},
    representations::Bispinor,
};
use spenso::{
    network::{
        MinResultRank, Network, Sequential,
        library::{
            LibraryTensor,
            function_lib::{INBUILTS, Panic},
            symbolic::{ExplicitKey, TensorLibrary},
        },
        parsing::{ParseSettings as TensorParseSettings, ShadowedStructure, ShorthandParsing},
        store::NetworkStore,
    },
    structure::{
        PermutedStructure,
        representation::{LibraryRep, Minkowski, RepName},
        slot::{AbsInd, DummyAind, ParseableAind, SlotError},
    },
    tensors::{
        data::{DataTensor, StorageTensor},
        parametric::ParamTensor,
    },
};
use spenso_hep_lib::{gamma_data_weyl, gamma5_weyl_data, proj_m_data_weyl, proj_p_data_weyl};
#[cfg(test)]
use symbolica::parser::ParseSettings;
use symbolica::{
    atom::{
        Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol, representation::FunView,
    },
    function,
    graph::Graph,
    id::Replacement,
    symbol,
};
use thiserror::Error;

use crate::{DiagramGroup, GraphGroupingOptions, GroupMember, NumeratorGrouping};

#[derive(Debug, Error)]
pub enum GroupingError {
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error("failed to construct a topology grouping key: {0}")]
    TopologyConstruction(String),
    #[error("source diagram index {0} does not fit in a symbolic integer")]
    SourceDiagramIndexOverflow(usize),
    #[error(
        "failed to merge physical cuts from group members {members:?} into master diagram '{master}': {source}"
    )]
    CutMerge {
        master: String,
        members: Vec<String>,
        #[source]
        source: Box<DiagramError>,
    },
    #[error(
        "failed to merge topology threshold candidates from group members {members:?} into master diagram '{master}': {source}"
    )]
    ThresholdCandidateMerge {
        master: String,
        members: Vec<String>,
        #[source]
        source: Box<DiagramError>,
    },
    #[error(
        "failed to transport grouped partition from member diagram '{member}' to master diagram '{master}': {message}"
    )]
    PartitionTransport {
        member: String,
        master: String,
        message: String,
    },
    #[error(
        "failed to evaluate numerator tensors for diagram '{diagram}', sample {sample}: {message}"
    )]
    TensorEvaluation {
        diagram: String,
        sample: usize,
        message: String,
    },
}

/// Typed abstract indices used while comparing finalized FeynKit numerators.
///
/// FeynKit deliberately keeps endpoint and owner identity in public
/// expressions.  Spenso's tensor canonicalizer needs to parse those atoms in
/// order to alpha-rename *contracted* indices while leaving dangling endpoint
/// indices untouched.  Keeping this type here makes numerator normalization a
/// FeynKit operation instead of reintroducing GammaLoop's historical `Aind`
/// adapter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum GroupingIndex {
    Normal(usize),
    Symbol(Symbol),
    Source(usize, usize),
    Sink(usize, usize),
    Vertex(usize, usize),
    Edge(usize, usize),
    Dummy(usize),
}

static GROUPING_DUMMY_COUNTER: AtomicUsize = AtomicUsize::new(0);

impl AbsInd for GroupingIndex {}

impl DummyAind for GroupingIndex {
    fn new_dummy() -> Self {
        Self::Dummy(GROUPING_DUMMY_COUNTER.fetch_add(1, Ordering::Relaxed))
    }

    fn new_dummy_at(index: usize) -> Self {
        Self::Dummy(index)
    }

    fn is_dummy(&self) -> bool {
        matches!(self, Self::Dummy(_))
    }
}

impl fmt::Display for GroupingIndex {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Normal(index) => index.fmt(formatter),
            Self::Symbol(symbol) => symbol.fmt(formatter),
            Self::Source(owner, local) => write!(formatter, "source({owner},{local})"),
            Self::Sink(owner, local) => write!(formatter, "sink({owner},{local})"),
            Self::Vertex(owner, local) => write!(formatter, "vertex({owner},{local})"),
            Self::Edge(owner, local) => write!(formatter, "edge({owner},{local})"),
            Self::Dummy(index) => write!(formatter, "dummy({index})"),
        }
    }
}

#[derive(Debug, Error)]
enum GroupingIndexError {
    #[error("index must be a non-negative integer, got {0}")]
    NotNatural(String),
    #[error("unsupported FeynKit index {0}")]
    Unsupported(String),
    #[error("{name} expects two non-negative integer arguments, got {actual}")]
    InvalidArity { name: &'static str, actual: usize },
}

impl From<GroupingIndexError> for SlotError {
    fn from(error: GroupingIndexError) -> Self {
        SlotError::Any(error.into())
    }
}

fn natural_index(view: AtomView<'_>) -> Result<usize, GroupingIndexError> {
    i64::try_from(view)
        .ok()
        .and_then(|value| usize::try_from(value).ok())
        .ok_or_else(|| GroupingIndexError::NotNatural(view.to_string()))
}

fn owned_index(
    function: FunView<'_>,
    name: &'static str,
    constructor: impl FnOnce(usize, usize) -> GroupingIndex,
) -> Result<GroupingIndex, GroupingIndexError> {
    if function.get_nargs() != 2 {
        return Err(GroupingIndexError::InvalidArity {
            name,
            actual: function.get_nargs(),
        });
    }
    let mut arguments = function.iter();
    let owner = natural_index(arguments.next().expect("arity checked"))?;
    let local = natural_index(arguments.next().expect("arity checked"))?;
    Ok(constructor(owner, local))
}

impl ParseableAind for GroupingIndex {
    type Error = GroupingIndexError;

    fn from_view(view: AtomView<'_>) -> Result<Self, Self::Error> {
        match view {
            AtomView::Num(_) => natural_index(view).map(Self::Normal),
            AtomView::Var(variable) => Ok(Self::Symbol(variable.get_symbol())),
            AtomView::Fun(function) => {
                let head = function.get_symbol();
                if head == symbol!("FeynKit::SourceIndex") {
                    owned_index(function, "FeynKit::SourceIndex", Self::Source)
                } else if head == symbol!("FeynKit::SinkIndex") {
                    owned_index(function, "FeynKit::SinkIndex", Self::Sink)
                } else if head == symbol!("FeynKit::VertexDummy") {
                    owned_index(function, "FeynKit::VertexDummy", Self::Vertex)
                } else if head == symbol!("FeynKit::EdgeDummy") {
                    owned_index(function, "FeynKit::EdgeDummy", Self::Edge)
                } else if head == symbol!("FeynKit::DummyIndex") {
                    if function.get_nargs() != 1 {
                        return Err(GroupingIndexError::InvalidArity {
                            name: "FeynKit::DummyIndex",
                            actual: function.get_nargs(),
                        });
                    }
                    let index = natural_index(function.iter().next().expect("arity checked"))?;
                    GROUPING_DUMMY_COUNTER.fetch_max(index + 1, Ordering::Relaxed);
                    Ok(Self::Dummy(index))
                } else {
                    Err(GroupingIndexError::Unsupported(view.to_string()))
                }
            }
            _ => Err(GroupingIndexError::Unsupported(view.to_string())),
        }
    }

    fn to_atom(&self) -> Atom {
        (*self).into()
    }
}

impl From<GroupingIndex> for Atom {
    fn from(index: GroupingIndex) -> Self {
        match index {
            GroupingIndex::Normal(index) => Atom::num(index as i64),
            GroupingIndex::Symbol(symbol) => Atom::var(symbol),
            GroupingIndex::Source(owner, local) => {
                function!(symbol!("FeynKit::SourceIndex"), owner as i64, local as i64)
            }
            GroupingIndex::Sink(owner, local) => {
                function!(symbol!("FeynKit::SinkIndex"), owner as i64, local as i64)
            }
            GroupingIndex::Vertex(owner, local) => {
                function!(symbol!("FeynKit::VertexDummy"), owner as i64, local as i64)
            }
            GroupingIndex::Edge(owner, local) => {
                function!(symbol!("FeynKit::EdgeDummy"), owner as i64, local as i64)
            }
            GroupingIndex::Dummy(index) => {
                function!(symbol!("FeynKit::DummyIndex"), index as i64)
            }
        }
    }
}

impl From<GroupingIndex> for AtomOrView<'_> {
    fn from(index: GroupingIndex) -> Self {
        Self::Atom(index.into())
    }
}

pub(crate) struct GroupingOutcome {
    pub diagrams: Vec<FeynmanDiagram>,
    pub groups: Vec<DiagramGroup>,
    pub zero_numerator_count: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum TopologyVertex {
    Internal,
    External(ExternalLeg),
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum InternalParticleColor {
    MassAndSpin { mass: ParameterId, spin: i64 },
    Exact(ParticleId),
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum TopologyEdge {
    Internal(InternalParticleColor),
    External {
        particle: ParticleId,
        directed: bool,
    },
}

type TopologyKey = Graph<TopologyVertex, TopologyEdge>;

struct PreparedNumerator {
    exact: Atom,
    sample_source: Atom,
    samples: Vec<Atom>,
}

type GroupingTensorLibrary = TensorLibrary<ParamTensor<ExplicitKey<GroupingIndex>>, GroupingIndex>;
type GroupingTensorNetwork = Network<
    NetworkStore<ParamTensor<ShadowedStructure<GroupingIndex>>, Atom>,
    ExplicitKey<GroupingIndex>,
    Symbol,
    GroupingIndex,
>;

#[derive(Clone)]
struct CanonicalTopology {
    key: TopologyKey,
    /// Mapping from this diagram's input vertex IDs into `key`'s canonical IDs.
    vertex_map: Vec<usize>,
    /// Whether the canonical key was obtained through the global incoming/outgoing swap.
    left_right_swapped: bool,
}

#[derive(Clone, Copy)]
enum ComparisonMode {
    Identical,
    UpToSign,
    UpToScalar,
}

pub(crate) fn group_diagrams(
    diagrams: Vec<FeynmanDiagram>,
    model: &Model,
    grouping: &NumeratorGrouping,
    symmetrize_left_right: bool,
) -> Result<GroupingOutcome, GroupingError> {
    if matches!(grouping, NumeratorGrouping::None) {
        let groups = singleton_groups(&diagrams, &(0..diagrams.len()).collect::<Vec<_>>());
        return Ok(GroupingOutcome {
            diagrams,
            groups,
            zero_numerator_count: 0,
        });
    }

    let scalar_names = scalar_names(model);
    let compares_numerators = matches!(
        grouping,
        NumeratorGrouping::Identical(_)
            | NumeratorGrouping::UpToSign(_)
            | NumeratorGrouping::UpToScalar(_)
    );
    let compares_canonical_numerators = match grouping {
        NumeratorGrouping::Identical(options)
        | NumeratorGrouping::UpToSign(options)
        | NumeratorGrouping::UpToScalar(options) => options.check_canonical_numerator,
        NumeratorGrouping::None | NumeratorGrouping::OnlyDetectZeroes => false,
    };
    let symmetric_polarizations = match grouping {
        NumeratorGrouping::Identical(options)
        | NumeratorGrouping::UpToSign(options)
        | NumeratorGrouping::UpToScalar(options) => options.symmetric_polarizations,
        NumeratorGrouping::None | NumeratorGrouping::OnlyDetectZeroes => false,
    };
    let mut retained = Vec::with_capacity(diagrams.len());
    let mut source_diagrams = Vec::with_capacity(diagrams.len());
    let mut prepared = Vec::with_capacity(diagrams.len());
    let mut zero_numerator_count = 0;
    for (source_diagram, diagram) in diagrams.into_iter().enumerate() {
        // Contract color before abstract-index canonicalization. Canonizing a
        // raw product of color tensors can encounter the same concrete base
        // slot more than once, while the projector closes precisely those
        // external slots.
        let complete = model
            .expand_couplings(
                &(diagram.numerator() * diagram.numerator_prefactor() * diagram.projector()),
            )
            .to_parametric_color()
            .simplify_color_with(ColorSimplifySettings::default().with_cof_dimension_invariants());
        let zero_check = complete
            .expand_color()
            .into_iter()
            .fold(Atom::Zero, |sum, (color, lorentz)| sum + color * lorentz);
        if zero_check.expand().is_zero() {
            zero_numerator_count += 1;
            continue;
        }
        let sample_source = if compares_numerators {
            // Preserve the symbolically simplified color tensor here.  Color
            // expansion is useful for the zero proof above and for tensor
            // samples, but it resolves representation dimensions to their
            // concrete values.  Feeding that expanded form to the exact path
            // would turn invariants such as `Nc` into `3` before an exact
            // rescaling ratio can be extracted.
            normalize_momentum_routing(&diagram, complete)
        } else {
            Atom::one()
        };
        // Exact symbolic comparison needs alpha-renamed dummy indices.  Keep
        // the routed, uncanonized expression separately for numerical tensor
        // execution: canonicalization is an optional comparison strategy and
        // historically was not applied to the numerical sample network.
        let exact = if compares_canonical_numerators {
            let exact_source = if symmetric_polarizations {
                normalize_symmetric_polarizations(&sample_source)
            } else {
                sample_source.clone()
            };
            exact_source.canonize(GroupingIndex::Dummy)
        } else {
            Atom::one()
        };
        retained.push(diagram);
        source_diagrams.push(source_diagram);
        prepared.push(PreparedNumerator {
            exact,
            sample_source,
            samples: Vec::new(),
        });
    }

    let (options, mode) = match grouping {
        NumeratorGrouping::Identical(options) => (options, ComparisonMode::Identical),
        NumeratorGrouping::UpToSign(options) => (options, ComparisonMode::UpToSign),
        NumeratorGrouping::UpToScalar(options) => (options, ComparisonMode::UpToScalar),
        NumeratorGrouping::None | NumeratorGrouping::OnlyDetectZeroes => {
            return Ok(GroupingOutcome {
                groups: singleton_groups(&retained, &source_diagrams),
                diagrams: retained,
                zero_numerator_count,
            });
        }
    };
    if matches!(mode, ComparisonMode::UpToScalar) {
        for numerator in &mut prepared {
            numerator.exact = numerator.exact.collect_factors();
        }
    }

    let mut buckets = BTreeMap::<TopologyKey, Vec<usize>>::new();
    let mut canonical_frames = Vec::with_capacity(retained.len());
    for (index, (diagram, numerator)) in retained.iter().zip(&mut prepared).enumerate() {
        let canonical = canonical_topology(diagram, model, options, symmetrize_left_right)?;
        let key = canonical.key.clone();
        canonical_frames.push(canonical);
        numerator.samples = numerical_tensor_samples(&numerator.sample_source, diagram, options)?;
        if !numerator.samples.is_empty()
            && numerator
                .samples
                .iter()
                .all(|sample| sample.expand().is_zero())
        {
            zero_numerator_count += 1;
            continue;
        }
        buckets.entry(key).or_default().push(index);
    }
    let mut groups = Vec::new();
    for indices in buckets.into_values() {
        let mut bucket_groups: Vec<DiagramGroup> = Vec::new();
        for diagram in indices {
            let matching_group =
                bucket_groups
                    .iter()
                    .enumerate()
                    .find_map(|(group_index, group)| {
                        compare(
                            &prepared[diagram],
                            &prepared[group.master],
                            mode,
                            options,
                            &scalar_names,
                        )
                        .map(|ratio| (group_index, ratio))
                    });
            if let Some((group_index, ratio)) = matching_group {
                bucket_groups[group_index].members.push(GroupMember {
                    source_diagram: source_diagrams[diagram],
                    source_id: retained[diagram].id(),
                    source_name: retained[diagram].name().to_owned(),
                    diagram,
                    ratio,
                    overall_factor: retained[diagram].overall_factor().clone(),
                });
            } else {
                bucket_groups.push(DiagramGroup {
                    master: diagram,
                    members: vec![GroupMember {
                        source_diagram: source_diagrams[diagram],
                        source_id: retained[diagram].id(),
                        source_name: retained[diagram].name().to_owned(),
                        diagram,
                        ratio: Atom::one(),
                        overall_factor: retained[diagram].overall_factor().clone(),
                    }],
                });
            }
        }
        groups.extend(bucket_groups);
    }
    groups.sort_by_key(|group| group.master);

    let (diagrams, groups) = collapse_groups(retained, &canonical_frames, groups)?;
    Ok(GroupingOutcome {
        diagrams,
        groups,
        zero_numerator_count,
    })
}

/// Express every edge momentum in a diagram-independent loop/external basis.
///
/// Public diagram numerators intentionally use stable edge IDs. Those IDs are
/// ideal for inspecting one diagram, but isomorphic diagrams can route the
/// same momentum through different edges. Grouping instead labels independent
/// loop momenta by basis position and external momenta by physical connection,
/// matching the normalization that historically preceded numerator
/// comparison.
fn normalize_momentum_routing(diagram: &FeynmanDiagram, mut numerator: Atom) -> Atom {
    let basis = diagram.loop_momentum_basis();
    let external_connections = diagram
        .edges()
        .filter_map(|(edge, endpoints, _)| {
            diagram
                .vertex(endpoints.source)
                .and_then(|vertex| vertex.external.as_ref())
                .or_else(|| {
                    diagram
                        .vertex(endpoints.target)
                        .and_then(|vertex| vertex.external.as_ref())
                })
                .map(|external| (edge, external.connection))
        })
        .collect::<BTreeMap<_, _>>();
    let polarization_symbols = [PS.u, PS.ubar, PS.v, PS.vbar, PS.eps, PS.ebar];

    numerator = numerator.replace_map(|term, _, output| {
        let AtomView::Fun(function) = term else {
            return;
        };
        if function.get_symbol() == crate::momentum_symbol()
            && (1..=2).contains(&function.get_nargs())
        {
            let mut arguments = function.iter();
            let Some(edge) = arguments.next().and_then(momentum_edge) else {
                return;
            };
            let Some(signature) = basis.edge_signatures.get(&edge) else {
                return;
            };
            let tensor_index = arguments.next().map(|argument| argument.to_owned());
            let mut replacement = Atom::Zero;
            for (loop_index, coefficient) in signature
                .loops
                .integer_coefficients()
                .into_iter()
                .enumerate()
            {
                if coefficient != 0 {
                    replacement += coefficient
                        * routed_momentum(
                            symbol!("FeynKit::LoopMomentum").call(loop_index as i64),
                            tensor_index.as_ref(),
                        );
                }
            }
            for (external_position, coefficient) in signature
                .external
                .integer_coefficients()
                .into_iter()
                .enumerate()
            {
                if coefficient == 0 {
                    continue;
                }
                let Some(external_edge) = basis.external_edges.get(external_position) else {
                    return;
                };
                let Some(connection) = external_connections.get(external_edge) else {
                    return;
                };
                replacement += coefficient
                    * routed_momentum(
                        symbol!("FeynKit::ExternalMomentum").call(*connection as i64),
                        tensor_index.as_ref(),
                    );
            }
            **output = replacement;
        } else if polarization_symbols.contains(&function.get_symbol()) && function.get_nargs() == 2
        {
            let mut arguments = function.iter();
            let Some(edge) = arguments.next().and_then(momentum_edge) else {
                return;
            };
            let Some(connection) = external_connections.get(&edge) else {
                return;
            };
            let index = arguments.next().expect("polarization arity checked");
            **output = FunctionBuilder::new(function.get_symbol())
                .add_arg(symbol!("FeynKit::ExternalMomentum").call(*connection as i64))
                .add_arg(index)
                .finish();
        }
    });
    numerator
}

/// Give left/right partners the same symbolic polarization tensor for exact
/// numerator comparisons, matching the shared numerical samples selected by
/// `symmetric_polarizations`.
fn normalize_symmetric_polarizations(numerator: &Atom) -> Atom {
    numerator.replace_map(|term, _, output| {
        let AtomView::Fun(function) = term else {
            return;
        };
        let normalized_head = if [PS.eps, PS.ebar].contains(&function.get_symbol()) {
            PS.eps
        } else if [PS.u, PS.ubar, PS.v, PS.vbar].contains(&function.get_symbol()) {
            PS.u
        } else {
            return;
        };
        **output = FunctionBuilder::new(normalized_head)
            .add_args(function.iter())
            .finish();
    })
}

fn momentum_edge(argument: AtomView<'_>) -> Option<EdgeId> {
    i64::try_from(argument)
        .ok()
        .and_then(|edge| usize::try_from(edge).ok())
        .map(EdgeId)
}

fn routed_momentum(label: Atom, tensor_index: Option<&Atom>) -> Atom {
    let mut momentum = FunctionBuilder::new(crate::momentum_symbol()).add_arg(label);
    if let Some(index) = tensor_index {
        momentum = momentum.add_arg(index);
    }
    momentum.finish()
}

fn partition_vertices(
    diagram: &FeynmanDiagram,
    half_edges: &[DiagramHalfEdge],
    master: &FeynmanDiagram,
) -> Result<BTreeSet<VertexId>, GroupingError> {
    let endpoints = diagram
        .edges()
        .map(|(edge, endpoints, _)| (edge, endpoints))
        .collect::<BTreeMap<_, _>>();
    half_edges
        .iter()
        .map(|half_edge| {
            let endpoints = endpoints.get(&half_edge.edge).ok_or_else(|| {
                GroupingError::PartitionTransport {
                    member: diagram.name().to_owned(),
                    master: master.name().to_owned(),
                    message: format!("unknown source edge {}", half_edge.edge.0),
                }
            })?;
            Ok(match half_edge.endpoint {
                DiagramEndpoint::Source => endpoints.source,
                DiagramEndpoint::Target => endpoints.target,
            })
        })
        .collect()
}

fn remap_partition(
    member: &FeynmanDiagram,
    member_frame: &CanonicalTopology,
    master: &FeynmanDiagram,
    master_frame: &CanonicalTopology,
    half_edges: &[DiagramHalfEdge],
) -> Result<Vec<DiagramHalfEdge>, GroupingError> {
    if member_frame.vertex_map.len() != member.vertices().count()
        || master_frame.vertex_map.len() != master.vertices().count()
    {
        return Err(GroupingError::PartitionTransport {
            member: member.name().to_owned(),
            master: master.name().to_owned(),
            message: "canonical vertex map does not cover the diagram".to_owned(),
        });
    }

    let mut canonical_to_master = BTreeMap::new();
    for (vertex, canonical) in master_frame.vertex_map.iter().copied().enumerate() {
        if canonical_to_master
            .insert(canonical, VertexId(vertex))
            .is_some()
        {
            return Err(GroupingError::PartitionTransport {
                member: member.name().to_owned(),
                master: master.name().to_owned(),
                message: format!("canonical master vertex {canonical} is not unique"),
            });
        }
    }

    let source_vertices = partition_vertices(member, half_edges, master)?;
    let mapped_vertices = source_vertices
        .into_iter()
        .map(|vertex| {
            let canonical = member_frame.vertex_map.get(vertex.0).ok_or_else(|| {
                GroupingError::PartitionTransport {
                    member: member.name().to_owned(),
                    master: master.name().to_owned(),
                    message: format!("source vertex {} has no canonical image", vertex.0),
                }
            })?;
            canonical_to_master.get(canonical).copied().ok_or_else(|| {
                GroupingError::PartitionTransport {
                    member: member.name().to_owned(),
                    master: master.name().to_owned(),
                    message: format!("canonical vertex {canonical} has no master preimage"),
                }
            })
        })
        .collect::<Result<BTreeSet<_>, _>>()?;

    Ok(master
        .edges()
        .flat_map(|(edge, endpoints, _)| {
            [
                mapped_vertices
                    .contains(&endpoints.source)
                    .then_some(DiagramHalfEdge {
                        edge,
                        endpoint: DiagramEndpoint::Source,
                    }),
                mapped_vertices
                    .contains(&endpoints.target)
                    .then_some(DiagramHalfEdge {
                        edge,
                        endpoint: DiagramEndpoint::Target,
                    }),
            ]
        })
        .flatten()
        .collect())
}

fn remap_partitions(
    member: &FeynmanDiagram,
    member_frame: &CanonicalTopology,
    master: &FeynmanDiagram,
    master_frame: &CanonicalTopology,
    left: &[DiagramHalfEdge],
    right: &[DiagramHalfEdge],
) -> Result<(Vec<DiagramHalfEdge>, Vec<DiagramHalfEdge>), GroupingError> {
    let mut left = remap_partition(member, member_frame, master, master_frame, left)?;
    let mut right = remap_partition(member, member_frame, master, master_frame, right)?;
    if member_frame.left_right_swapped != master_frame.left_right_swapped {
        std::mem::swap(&mut left, &mut right);
    }
    Ok((left, right))
}

fn collapse_groups(
    diagrams: Vec<FeynmanDiagram>,
    canonical_frames: &[CanonicalTopology],
    groups: Vec<DiagramGroup>,
) -> Result<(Vec<FeynmanDiagram>, Vec<DiagramGroup>), GroupingError> {
    let mut masters = Vec::with_capacity(groups.len());
    let mut collapsed_groups = Vec::with_capacity(groups.len());
    for group in groups {
        let output_index = masters.len();
        let topology_master = &diagrams[group.master];
        let master_frame = &canonical_frames[group.master];
        let master = topology_master.clone();
        // Numerator grouping is intentionally based on the denominator
        // topology. Stable half-edge IDs are local to each input diagram, so
        // transport the member's vertex partition through the canonical
        // topology frame before rebuilding cut endpoints on the master.
        let mut cut_partitions = Vec::new();
        for member in &group.members {
            let source = &diagrams[member.diagram];
            let source_frame = &canonical_frames[member.diagram];
            for cut in source.cuts() {
                cut_partitions.push(remap_partitions(
                    source,
                    source_frame,
                    topology_master,
                    master_frame,
                    &cut.left.half_edges,
                    &cut.right.half_edges,
                )?);
            }
        }
        cut_partitions.sort();
        cut_partitions.dedup();
        let master = if !cut_partitions.is_empty() || !master.cuts().is_empty() {
            let members = group
                .members
                .iter()
                .map(|member| diagrams[member.diagram].name().to_owned())
                .collect::<Vec<_>>();
            let master_name = master.name().to_owned();
            master
                .with_cut_partitions(cut_partitions)
                .map_err(|source| GroupingError::CutMerge {
                    master: master_name,
                    members,
                    source: Box::new(source),
                })?
        } else {
            master
        };
        let mut threshold_partitions = Vec::new();
        for member in &group.members {
            let source = &diagrams[member.diagram];
            let source_frame = &canonical_frames[member.diagram];
            for candidate in source.topology_threshold_candidates() {
                threshold_partitions.push(remap_partitions(
                    source,
                    source_frame,
                    topology_master,
                    master_frame,
                    &candidate.left,
                    &candidate.right,
                )?);
            }
        }
        threshold_partitions.sort();
        threshold_partitions.dedup();
        let master = if !threshold_partitions.is_empty()
            || !master.topology_threshold_candidates().is_empty()
        {
            let members = group
                .members
                .iter()
                .map(|member| diagrams[member.diagram].name().to_owned())
                .collect::<Vec<_>>();
            let master_name = master.name().to_owned();
            master
                .with_topology_threshold_partitions(threshold_partitions)
                .map_err(|source| GroupingError::ThresholdCandidateMerge {
                    master: master_name,
                    members,
                    source: Box::new(source),
                })?
        } else {
            master
        };
        let factor = if group.members.len() == 1 {
            master.overall_factor().clone()
        } else {
            let mut factor = Atom::Zero;
            for member in &group.members {
                let source_diagram = i64::try_from(member.source_diagram).map_err(|_| {
                    GroupingError::SourceDiagramIndexOverflow(member.source_diagram)
                })?;
                factor += function!(
                    symbol!("feynkit_generator_factor::NumeratorDependentGrouping"),
                    Atom::num(source_diagram),
                    member.ratio.clone(),
                    member.overall_factor.clone()
                );
            }
            factor
        };
        masters.push(master.with_overall_factor(factor));
        collapsed_groups.push(DiagramGroup {
            master: output_index,
            members: group
                .members
                .into_iter()
                .map(|mut member| {
                    member.diagram = output_index;
                    member
                })
                .collect(),
        });
    }
    Ok((masters, collapsed_groups))
}

fn singleton_groups(diagrams: &[FeynmanDiagram], source_diagrams: &[usize]) -> Vec<DiagramGroup> {
    diagrams
        .iter()
        .zip(source_diagrams.iter().copied())
        .enumerate()
        .map(|(diagram, (source, source_diagram))| DiagramGroup {
            master: diagram,
            members: vec![GroupMember {
                source_diagram,
                source_id: source.id(),
                source_name: source.name().to_owned(),
                diagram,
                ratio: Atom::one(),
                overall_factor: source.overall_factor().clone(),
            }],
        })
        .collect()
}

fn scalar_names(model: &Model) -> BTreeSet<String> {
    model
        .parameters()
        .iter()
        .map(|parameter| parameter.name.clone())
        .chain(
            model
                .couplings()
                .iter()
                .map(|coupling| coupling.name.clone()),
        )
        .chain(
            model
                .functions()
                .iter()
                .map(|function| function.name.clone()),
        )
        .chain(
            model
                .form_factors()
                .iter()
                .map(|form_factor| form_factor.name.clone()),
        )
        .collect()
}

fn numerical_tensor_samples(
    numerator: &Atom,
    diagram: &FeynmanDiagram,
    options: &GraphGroupingOptions,
) -> Result<Vec<Atom>, GroupingError> {
    (0..options.number_of_numerical_samples)
        .map(|sample| evaluate_tensor_sample(numerator, diagram, options, sample))
        .collect()
}

fn evaluate_tensor_sample(
    numerator: &Atom,
    diagram: &FeynmanDiagram,
    options: &GraphGroupingOptions,
    sample: usize,
) -> Result<Atom, GroupingError> {
    let library = tensor_sample_library(diagram, options, sample)?;
    let settings = TensorParseSettings {
        shorthand_parsing: ShorthandParsing::expand_schoonschip_only(),
        ..TensorParseSettings::default()
    };
    let scalar_replacements = if options.fully_numerical_substitution {
        diagram
            .model()
            .parameters()
            .iter()
            .filter(|parameter| parameter.value.is_some())
            .map(|parameter| {
                Replacement::new(
                    Atom::var(symbol!(&format!("UFO::{}", parameter.name))),
                    complex_sample_value(
                        &format!("model:UFO::{}", parameter.name),
                        options.numerical_sample_seed,
                        sample,
                    ),
                )
            })
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };

    let mut functions = Panic::new_lib();
    functions.insert(
        INBUILTS.conj,
        |tensor: ParamTensor<ShadowedStructure<GroupingIndex>>| {
            tensor.map_data_self(|value| value.conj())
        },
    );
    let color_replacements = [
        Replacement::new(Atom::var(CS.nc), Atom::num(3)),
        Replacement::new(Atom::var(CS.tr), Atom::num(1) / Atom::num(2)),
    ];
    numerator
        .expand_color()
        .into_iter()
        .try_fold(Atom::Zero, |sum, (color, lorentz)| {
            let mut network =
                GroupingTensorNetwork::try_from_view(lorentz.as_view(), &library, &settings)
                    .map_err(|error| tensor_evaluation_error(diagram, sample, error))?;
            for scalar in &mut network.store.scalar {
                *scalar = scalar.replace_multiple(&scalar_replacements);
            }
            network
                .execute::<Sequential, MinResultRank, _, _, _>(&library, &functions)
                .map_err(|error| tensor_evaluation_error(diagram, sample, error))?;
            let scalar: Atom = network
                .result_scalar()
                .map_err(|error| tensor_evaluation_error(diagram, sample, error))?
                .into();
            let color = color.canonize(GroupingIndex::Dummy);
            Ok(sum + (color * scalar).replace_multiple(&color_replacements))
        })
        .map(|evaluated| evaluated.expand())
}

fn tensor_sample_library(
    diagram: &FeynmanDiagram,
    options: &GraphGroupingOptions,
    sample: usize,
) -> Result<GroupingTensorLibrary, GroupingError> {
    let mut library = GroupingTensorLibrary::new();
    library.update_ids();
    insert_dirac_tensors(&mut library);

    for loop_index in 0..diagram.loop_momentum_basis().loop_edges.len() {
        insert_sample_vector(
            &mut library,
            crate::momentum_symbol(),
            symbol!("FeynKit::LoopMomentum").call(loop_index as i64),
            &format!("loop-momentum:{loop_index}"),
            false,
            options,
            sample,
        )?;
    }

    let mut external_legs = diagram
        .edges()
        .filter_map(|(_, endpoints, edge)| {
            diagram
                .vertex(endpoints.source)
                .and_then(|vertex| vertex.external.as_ref())
                .or_else(|| {
                    diagram
                        .vertex(endpoints.target)
                        .and_then(|vertex| vertex.external.as_ref())
                })
                .map(|external| (external.clone(), edge.particle))
        })
        .collect::<Vec<_>>();
    external_legs.sort_by_key(|(external, particle)| {
        (external.connection, external.state, particle.index())
    });

    let mut inserted_momenta = BTreeSet::new();
    for (external, particle_id) in external_legs {
        if inserted_momenta.insert(external.connection) {
            insert_sample_vector(
                &mut library,
                crate::momentum_symbol(),
                symbol!("FeynKit::ExternalMomentum").call(external.connection as i64),
                &format!("external-momentum:{}", external.connection),
                false,
                options,
                sample,
            )?;
        }

        let particle = diagram.model().particle_by_id(particle_id)?;
        let head = match (particle.spin, external.state, particle.is_antiparticle()) {
            (2, ExternalState::Incoming, false) => PS.u,
            (2, ExternalState::Incoming, true) => PS.vbar,
            (2, ExternalState::Outgoing, false) => PS.ubar,
            (2, ExternalState::Outgoing, true) => PS.v,
            (3, ExternalState::Incoming, _) => PS.eps,
            (3, ExternalState::Outgoing, _) => PS.ebar,
            _ => continue,
        };
        let label = if options.symmetric_polarizations {
            format!("polarization:{}", external.connection)
        } else {
            format!("polarization:{head}:{}", external.connection)
        };
        insert_sample_wavefunction(
            &mut library,
            head,
            symbol!("FeynKit::ExternalMomentum").call(external.connection as i64),
            &label,
            particle.spin,
            options,
            sample,
        )?;
    }
    Ok(library)
}

fn insert_dirac_tensors(library: &mut GroupingTensorLibrary) {
    let gamma = gamma_data_weyl(
        AGS.gamma_strct::<GroupingIndex>(4),
        Atom::num(1),
        Atom::num(0),
    )
    .map_data(|value| value.re + Atom::i() * value.im);
    library.insert_explicit(PermutedStructure::identity(ParamTensor::composite(
        DataTensor::Sparse(gamma),
    )));
    let gamma5 = gamma5_weyl_data(
        AGS.gamma5_strct::<GroupingIndex>(4),
        Atom::num(1),
        Atom::num(0),
    )
    .map_data(|value| value.re + Atom::i() * value.im);
    library.insert_explicit(PermutedStructure::identity(ParamTensor::composite(
        DataTensor::Sparse(gamma5),
    )));
    let projm = proj_m_data_weyl(
        AGS.projm_strct::<GroupingIndex>(4),
        Atom::num(1),
        Atom::num(0),
    )
    .map_data(|value| value.re + Atom::i() * value.im);
    library.insert_explicit(PermutedStructure::identity(ParamTensor::composite(
        DataTensor::Sparse(projm),
    )));
    let projp = proj_p_data_weyl(
        AGS.projp_strct::<GroupingIndex>(4),
        Atom::num(1),
        Atom::num(0),
    )
    .map_data(|value| value.re + Atom::i() * value.im);
    library.insert_explicit(PermutedStructure::identity(ParamTensor::composite(
        DataTensor::Sparse(projp),
    )));
}

#[allow(clippy::too_many_arguments)]
fn insert_sample_vector(
    library: &mut GroupingTensorLibrary,
    head: Symbol,
    argument: Atom,
    label: &str,
    complex: bool,
    options: &GraphGroupingOptions,
    sample: usize,
) -> Result<(), GroupingError> {
    let key = ExplicitKey::from_iter([Minkowski {}.new_rep(4)], head, Some(vec![argument]));
    let data = (0..4)
        .map(|component| {
            let component_label = format!("{label}:{component}");
            if complex {
                complex_sample_value(&component_label, options.numerical_sample_seed, sample)
            } else {
                Atom::num(sample_value(
                    &component_label,
                    options.numerical_sample_seed,
                    sample,
                ))
            }
        })
        .collect();
    let tensor = ParamTensor::from_dense(key.structure, data).map_err(|error| {
        GroupingError::TensorEvaluation {
            diagram: "sample library".to_owned(),
            sample,
            message: error.to_string(),
        }
    })?;
    library.insert_explicit(PermutedStructure::identity(tensor));
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn insert_sample_wavefunction(
    library: &mut GroupingTensorLibrary,
    head: Symbol,
    argument: Atom,
    label: &str,
    spin: i64,
    options: &GraphGroupingOptions,
    sample: usize,
) -> Result<(), GroupingError> {
    let representation = if spin == 2 {
        Bispinor {}.new_rep(4).cast::<LibraryRep>()
    } else {
        Minkowski {}.new_rep(4).cast::<LibraryRep>()
    };
    let key = ExplicitKey::from_iter([representation], head, Some(vec![argument]));
    let data = (0..4)
        .map(|component| {
            complex_sample_value(
                &format!("{label}:{component}"),
                options.numerical_sample_seed,
                sample,
            )
        })
        .collect();
    let tensor = ParamTensor::from_dense(key.structure, data).map_err(|error| {
        GroupingError::TensorEvaluation {
            diagram: "sample library".to_owned(),
            sample,
            message: error.to_string(),
        }
    })?;
    library.insert_explicit(PermutedStructure::identity(tensor));
    Ok(())
}

fn complex_sample_value(expression: &str, seed: u16, sample: usize) -> Atom {
    Atom::num(sample_value(expression, seed, sample))
        + Atom::i()
            * Atom::num(sample_value(
                &format!("{expression}:imaginary"),
                seed,
                sample,
            ))
}

fn tensor_evaluation_error(
    diagram: &FeynmanDiagram,
    sample: usize,
    error: impl fmt::Display,
) -> GroupingError {
    GroupingError::TensorEvaluation {
        diagram: diagram.name().to_owned(),
        sample,
        message: error.to_string(),
    }
}

fn sample_value(expression: &str, seed: u16, sample: usize) -> i64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64
        ^ u64::from(seed)
        ^ (sample as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    for byte in expression.bytes() {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    let magnitude = (hash % 29 + 2) as i64;
    if hash & (1 << 63) == 0 {
        magnitude
    } else {
        -magnitude
    }
}

fn compare(
    candidate: &PreparedNumerator,
    master: &PreparedNumerator,
    mode: ComparisonMode,
    options: &GraphGroupingOptions,
    scalar_names: &BTreeSet<String>,
) -> Option<Atom> {
    if options.check_canonical_numerator
        && let Some(ratio) = compare_atoms(&candidate.exact, &master.exact, mode, scalar_names)
    {
        return Some(ratio);
    }
    let mut ratios = candidate
        .samples
        .iter()
        .zip(&master.samples)
        .map(|(candidate, master)| compare_atoms(candidate, master, mode, scalar_names));
    let first = ratios.next()??;
    if ratios.all(|ratio| ratio.is_some_and(|ratio| expressions_equal(&ratio, &first))) {
        Some(first)
    } else {
        None
    }
}

fn compare_atoms(
    candidate: &Atom,
    master: &Atom,
    mode: ComparisonMode,
    scalar_names: &BTreeSet<String>,
) -> Option<Atom> {
    if expressions_equal(candidate, master) {
        return Some(Atom::num(1));
    }
    if matches!(mode, ComparisonMode::UpToSign) && (candidate + master).expand().is_zero() {
        return Some(Atom::num(-1));
    }
    if !matches!(mode, ComparisonMode::UpToScalar) || master.is_zero() {
        return None;
    }
    let ratio = (candidate / master).cancel();
    // A cancelled quotient is the exact symbolic rescaling.  Re-expanding
    // `ratio * master` is not a stronger check for tensor expressions: dummy
    // indices can be alpha-equivalent while carrying different canonical
    // labels in separately prepared diagrams.  Requiring literal equality of
    // that reconstruction rejects valid ratios.  Instead, as in the legacy
    // grouping implementation, reject the quotient whenever any tensor head
    // remains and otherwise accept the exact scalar.
    if expression_is_scalar(&ratio, scalar_names) {
        Some(ratio)
    } else {
        None
    }
}

fn expressions_equal(left: &Atom, right: &Atom) -> bool {
    (left - right).expand().is_zero()
}

fn expression_is_scalar(expression: &Atom, scalar_names: &BTreeSet<String>) -> bool {
    let color_scalars = [CS.na, CS.nc, CS.ca, CS.cf, CS.tr, CS.cas, CS.idx, CS.gram];
    expression.get_all_symbols(true).iter().all(|symbol| {
        symbol.is_scalar()
            || color_scalars.contains(symbol)
            || scalar_names.contains(symbol.get_stripped_name())
    })
}

fn canonical_topology(
    diagram: &FeynmanDiagram,
    model: &Model,
    options: &GraphGroupingOptions,
    symmetrize_left_right: bool,
) -> Result<CanonicalTopology, GroupingError> {
    canonical_topologies(diagram, model, options, symmetrize_left_right)?
        .into_iter()
        .min_by(|left, right| {
            left.key
                .cmp(&right.key)
                .then(left.left_right_swapped.cmp(&right.left_right_swapped))
        })
        .ok_or_else(|| GroupingError::TopologyConstruction("no canonical topology".to_owned()))
}

fn canonical_topologies(
    diagram: &FeynmanDiagram,
    model: &Model,
    options: &GraphGroupingOptions,
    symmetrize_left_right: bool,
) -> Result<Vec<CanonicalTopology>, GroupingError> {
    let mut topologies = vec![topology_key_variant(diagram, model, options, None)?];
    if symmetrize_left_right && let Some(partners) = left_right_partners(diagram) {
        topologies.push(topology_key_variant(
            diagram,
            model,
            options,
            Some(&partners),
        )?);
    }
    Ok(topologies)
}

fn left_right_partners(diagram: &FeynmanDiagram) -> Option<BTreeMap<ExternalLeg, ExternalLeg>> {
    let legs = diagram
        .vertices()
        .filter_map(|(_, vertex)| vertex.external.clone())
        .collect::<Vec<_>>();
    let by_connection = legs
        .iter()
        .cloned()
        .map(|leg| ((leg.connection, leg.state), leg))
        .collect::<BTreeMap<_, _>>();
    if by_connection.len() != legs.len() {
        return None;
    }
    legs.into_iter()
        .map(|leg| {
            let opposite = match leg.state {
                feynkit_graph::ExternalState::Incoming => feynkit_graph::ExternalState::Outgoing,
                feynkit_graph::ExternalState::Outgoing => feynkit_graph::ExternalState::Incoming,
            };
            by_connection
                .get(&(leg.connection, opposite))
                .cloned()
                .map(|partner| (leg, partner))
        })
        .collect()
}

fn topology_key_variant(
    diagram: &FeynmanDiagram,
    model: &Model,
    options: &GraphGroupingOptions,
    left_right_partners: Option<&BTreeMap<ExternalLeg, ExternalLeg>>,
) -> Result<CanonicalTopology, GroupingError> {
    let mut graph = Graph::new();
    for (_, vertex) in diagram.vertices() {
        let external = vertex.external.clone().map(|leg| {
            left_right_partners
                .and_then(|partners| partners.get(&leg))
                .cloned()
                .unwrap_or(leg)
        });
        graph.add_node(external.map_or(TopologyVertex::Internal, TopologyVertex::External));
    }
    for (_, endpoints, edge) in diagram.edges() {
        let external = diagram
            .vertex(endpoints.source)
            .is_some_and(|vertex| vertex.is_external())
            || diagram
                .vertex(endpoints.target)
                .is_some_and(|vertex| vertex.is_external());
        let (directed, color) = if external {
            // External edges retain both their exact particle species and
            // direction: exchanging different mass-degenerate external states
            // would change the process rather than merely its topology.
            (
                edge.directed,
                TopologyEdge::External {
                    particle: edge.particle,
                    directed: edge.directed,
                },
            )
        } else {
            let particle = model.particle_by_id(edge.particle)?;
            // Spin remains part of the reduced internal color so massless
            // fermions and vectors cannot be interchanged by canonicalization.
            let color = if options.differentiate_particle_masses_only {
                InternalParticleColor::MassAndSpin {
                    mass: particle.mass,
                    spin: particle.spin,
                }
            } else {
                InternalParticleColor::Exact(edge.particle)
            };
            (false, TopologyEdge::Internal(color))
        };
        let (source, target) = if left_right_partners.is_some() && external && directed {
            (endpoints.target.0, endpoints.source.0)
        } else {
            (endpoints.source.0, endpoints.target.0)
        };
        graph
            .add_edge(source, target, directed, color)
            .map_err(|error| GroupingError::TopologyConstruction(error.to_string()))?;
    }
    let canonical = graph.canonize();
    Ok(CanonicalTopology {
        key: canonical.graph,
        vertex_map: canonical.vertex_map,
        left_right_swapped: left_right_partners.is_some(),
    })
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use feynkit_graph::{
        DiagramCut, DiagramCutSide, DiagramEdge, DiagramEndpoint, DiagramHalfEdge,
        DiagramThresholdCandidate, DiagramVertex, ExternalState, VertexId,
    };
    use spenso::structure::representation::{Minkowski, RepName};

    use super::*;

    fn test_atom(expression: &str) -> Atom {
        Atom::parse(
            expression,
            "feynkit_grouping_test",
            ParseSettings::default(),
        )
        .unwrap()
    }

    fn model() -> Model {
        Model::from_json(
            r#"{
                "name":"grouping",
                "restriction":null,
                "orders":[],
                "parameters":[
                    {"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null},
                    {"name":"M","lhablock":null,"lhacode":null,"nature":"external","parameter_type":"real","value":[1.0,0.0],"expression":null},
                    {"name":"N","lhablock":null,"lhacode":null,"nature":"external","parameter_type":"real","value":[2.0,0.0],"expression":null}
                ],
                "particles":[
                    {"pdg_code":25,"name":"phi","antiname":"phi","spin":1,"color":1,"mass":"M","width":"ZERO","texname":"phi","antitexname":"phi","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":35,"name":"chi","antiname":"chi","spin":1,"color":1,"mass":"M","width":"ZERO","texname":"chi","antitexname":"chi","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":45,"name":"eta","antiname":"eta","spin":1,"color":1,"mass":"N","width":"ZERO","texname":"eta","antitexname":"eta","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":21,"name":"A","antiname":"A","spin":3,"color":1,"mass":"M","width":"ZERO","texname":"A","antitexname":"A","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0}
                ],
                "propagators":[],
                "lorentz_structures":[
                    {"name":"L","spins":[1,1],"structure":"1"}
                ],
                "couplings":[
                    {"name":"GC","expression":"1","orders":[],"value":null}
                ],
                "vertex_rules":[
                    {"name":"V","particles":["phi","phi"],"color_structures":["1"],"lorentz_structures":["L"],"couplings":[["GC"]]}
                ],
                "functions":[],
                "form_factors":[]
            }"#,
        )
        .unwrap()
    }

    fn particle(model: &Model, pdg: i64) -> ParticleId {
        model.particle_id_by_pdg(pdg).unwrap()
    }

    fn diagram(
        model: &Arc<Model>,
        name: &str,
        numerator: &str,
        internal_pdg: i64,
        external_pdg: i64,
    ) -> FeynmanDiagram {
        diagram_with_projector(model, name, numerator, "1", internal_pdg, external_pdg)
    }

    fn diagram_with_projector(
        model: &Arc<Model>,
        name: &str,
        numerator: &str,
        projector: &str,
        internal_pdg: i64,
        external_pdg: i64,
    ) -> FeynmanDiagram {
        let numerator =
            Atom::parse(numerator, "feynkit_grouping_test", ParseSettings::default()).unwrap();
        let projector =
            Atom::parse(projector, "feynkit_grouping_test", ParseSettings::default()).unwrap();
        diagram_with_atoms(
            model,
            name,
            numerator,
            projector,
            internal_pdg,
            external_pdg,
        )
    }

    fn diagram_with_atoms(
        model: &Arc<Model>,
        name: &str,
        numerator: Atom,
        projector: Atom,
        internal_pdg: i64,
        external_pdg: i64,
    ) -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), name)
            .numerator(numerator)
            .projector(projector);
        builder.add_vertex(DiagramVertex::external("in", 0, ExternalState::Incoming));
        let interaction = model.vertex_rule_id("V").unwrap();
        builder.add_vertex(DiagramVertex::interaction("left", interaction));
        builder.add_vertex(DiagramVertex::interaction("right", interaction));
        builder.add_vertex(DiagramVertex::external("out", 1, ExternalState::Outgoing));
        builder
            .add_edge(
                VertexId(0),
                VertexId(1),
                DiagramEdge::new(particle(model, external_pdg), false),
            )
            .unwrap();
        builder
            .add_edge(
                VertexId(1),
                VertexId(2),
                DiagramEdge::new(particle(model, internal_pdg), false),
            )
            .unwrap();
        builder
            .add_edge(
                VertexId(2),
                VertexId(3),
                DiagramEdge::new(particle(model, external_pdg), false),
            )
            .unwrap();
        builder.build().unwrap()
    }

    fn routed_projector_diagram(
        model: &Arc<Model>,
        name: &str,
        internal_edges_first: bool,
    ) -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), name);
        let incoming = builder.add_vertex(DiagramVertex::external_in_connection(
            "in",
            0,
            ExternalState::Incoming,
            0,
        ));
        let interaction = model.vertex_rule_id("V").unwrap();
        let left = builder.add_vertex(DiagramVertex::interaction("left", interaction));
        let right = builder.add_vertex(DiagramVertex::interaction("right", interaction));
        let outgoing = builder.add_vertex(DiagramVertex::external_in_connection(
            "out",
            1,
            ExternalState::Outgoing,
            1,
        ));
        let scalar = || DiagramEdge::new(particle(model, 25), false);
        let vector = || DiagramEdge::new(particle(model, 21), false);
        let add_internals = |builder: &mut feynkit_graph::FeynmanDiagramBuilder| {
            builder.add_edge(left, right, scalar()).unwrap();
            builder.add_edge(left, right, scalar()).unwrap()
        };
        let (incoming_edge, loop_edge) = if internal_edges_first {
            let loop_edge = add_internals(&mut builder);
            let incoming_edge = builder.add_edge(incoming, left, vector()).unwrap();
            builder.add_edge(right, outgoing, vector()).unwrap();
            (incoming_edge, loop_edge)
        } else {
            let incoming_edge = builder.add_edge(incoming, left, vector()).unwrap();
            let loop_edge = add_internals(&mut builder);
            builder.add_edge(right, outgoing, vector()).unwrap();
            (incoming_edge, loop_edge)
        };
        let endpoint = symbol!("FeynKit::SinkIndex").call((incoming_edge.0 as i64, 1_i64));
        let index = Minkowski {}.new_rep(4).to_symbolic([endpoint]);
        let numerator = FunctionBuilder::new(crate::momentum_symbol())
            .add_arg(loop_edge.0 as i64)
            .add_arg(index.clone())
            .finish();
        let projector = FunctionBuilder::new(PS.eps)
            .add_arg(incoming_edge.0 as i64)
            .add_arg(index)
            .finish();
        let diagram = builder
            .numerator(numerator)
            .projector(projector)
            .build()
            .unwrap();
        assert_eq!(diagram.loop_momentum_basis().loop_edges, vec![loop_edge]);
        diagram
    }

    fn cut_diagram(model: &Arc<Model>, mirrored: bool) -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), "cut-line");
        let incoming = builder.add_vertex(DiagramVertex::external_in_connection(
            "in",
            0,
            ExternalState::Incoming,
            0,
        ));
        let interaction = builder.add_vertex(DiagramVertex::interaction(
            "interaction",
            model.vertex_rule_id("V").unwrap(),
        ));
        let outgoing = builder.add_vertex(DiagramVertex::external_in_connection(
            "out",
            1,
            ExternalState::Outgoing,
            0,
        ));
        let scalar = || DiagramEdge::new(particle(model, 25), false);
        let incoming_edge = builder.add_edge(incoming, interaction, scalar()).unwrap();
        let outgoing_edge = builder.add_edge(interaction, outgoing, scalar()).unwrap();
        let incoming_source = DiagramHalfEdge {
            edge: incoming_edge,
            endpoint: DiagramEndpoint::Source,
        };
        let incoming_target = DiagramHalfEdge {
            edge: incoming_edge,
            endpoint: DiagramEndpoint::Target,
        };
        let outgoing_source = DiagramHalfEdge {
            edge: outgoing_edge,
            endpoint: DiagramEndpoint::Source,
        };
        let outgoing_target = DiagramHalfEdge {
            edge: outgoing_edge,
            endpoint: DiagramEndpoint::Target,
        };
        let side = |half_edges| DiagramCutSide {
            half_edges,
            coupling_orders: BTreeMap::new(),
            loop_count: 0,
        };
        let mut cut = DiagramCut {
            cut: vec![incoming_source],
            left: side(vec![incoming_source]),
            right: side(vec![incoming_target, outgoing_source, outgoing_target]),
        };
        if mirrored {
            // The global transformation exchanges both physical sides at
            // once. This is intentionally not just an independent relabeling
            // of the cut nodes.
            cut.cut = vec![outgoing_source];
            cut.left = side(vec![incoming_source, incoming_target, outgoing_source]);
            cut.right = side(vec![outgoing_target]);
        }
        let diagram = builder.cuts(vec![cut]).build().unwrap();
        diagram.validate().unwrap();
        diagram
    }

    fn partitioned_diagram(
        model: &Arc<Model>,
        name: &str,
        permuted_ids: bool,
        global_left_right_swap: bool,
    ) -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), name);
        let internal = |name: &str| DiagramVertex {
            name: name.to_owned(),
            interaction: None,
            external: None,
            numerator: Atom::one(),
        };
        let (incoming, left, right, outgoing, marker) = if permuted_ids {
            let marker = builder.add_vertex(internal("marker"));
            let outgoing = builder.add_vertex(DiagramVertex::external_in_connection(
                "out",
                1,
                ExternalState::Outgoing,
                0,
            ));
            let right = builder.add_vertex(internal("right"));
            let incoming = builder.add_vertex(DiagramVertex::external_in_connection(
                "in",
                0,
                ExternalState::Incoming,
                0,
            ));
            let left = builder.add_vertex(internal("left"));
            (incoming, left, right, outgoing, marker)
        } else {
            let incoming = builder.add_vertex(DiagramVertex::external_in_connection(
                "in",
                0,
                ExternalState::Incoming,
                0,
            ));
            let left = builder.add_vertex(internal("left"));
            let right = builder.add_vertex(internal("right"));
            let outgoing = builder.add_vertex(DiagramVertex::external_in_connection(
                "out",
                1,
                ExternalState::Outgoing,
                0,
            ));
            let marker = builder.add_vertex(internal("marker"));
            (incoming, left, right, outgoing, marker)
        };
        let (source_side, target_side) = if global_left_right_swap {
            (right, left)
        } else {
            (left, right)
        };
        let scalar = || DiagramEdge::new(particle(model, 25), false);
        let mut edges = Vec::new();
        let mut add_edge = |source, target| {
            let edge = builder.add_edge(source, target, scalar()).unwrap();
            edges.push((edge, source, target));
        };
        if permuted_ids {
            add_edge(left, marker);
            add_edge(target_side, outgoing);
            add_edge(left, right);
            add_edge(incoming, source_side);
            add_edge(left, right);
        } else {
            add_edge(incoming, source_side);
            add_edge(left, right);
            add_edge(left, right);
            add_edge(left, marker);
            add_edge(target_side, outgoing);
        }

        let left_vertices = if global_left_right_swap {
            BTreeSet::from([incoming, right])
        } else {
            BTreeSet::from([incoming, left, marker])
        };
        let mut left_half_edges = Vec::new();
        let mut right_half_edges = Vec::new();
        let mut crossing = Vec::new();
        for (edge, source, target) in edges {
            let source_half = DiagramHalfEdge {
                edge,
                endpoint: DiagramEndpoint::Source,
            };
            let target_half = DiagramHalfEdge {
                edge,
                endpoint: DiagramEndpoint::Target,
            };
            let source_is_left = left_vertices.contains(&source);
            let target_is_left = left_vertices.contains(&target);
            if source_is_left {
                left_half_edges.push(source_half);
            } else {
                right_half_edges.push(source_half);
            }
            if target_is_left {
                left_half_edges.push(target_half);
            } else {
                right_half_edges.push(target_half);
            }
            if source_is_left != target_is_left {
                crossing.push(if source_is_left {
                    source_half
                } else {
                    target_half
                });
            }
        }
        let side = |half_edges| DiagramCutSide {
            half_edges,
            coupling_orders: BTreeMap::new(),
            loop_count: 0,
        };
        let cut = DiagramCut {
            cut: crossing.clone(),
            left: side(left_half_edges.clone()),
            right: side(right_half_edges.clone()),
        };
        let candidate = DiagramThresholdCandidate {
            cut: crossing,
            left: left_half_edges,
            right: right_half_edges,
        };
        let diagram = builder
            .cuts(vec![cut])
            .topology_threshold_candidates(vec![candidate])
            .build()
            .unwrap();
        diagram.validate().unwrap();
        diagram
    }

    #[derive(Clone, Copy)]
    enum LeftRightAttachment {
        Identity,
        GlobalSwap,
        PartialSwap,
    }

    fn left_right_diagram(model: &Arc<Model>, attachment: LeftRightAttachment) -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), "left-right");
        let external = [
            builder.add_vertex(DiagramVertex::external_in_connection(
                "in0",
                0,
                ExternalState::Incoming,
                0,
            )),
            builder.add_vertex(DiagramVertex::external_in_connection(
                "in1",
                1,
                ExternalState::Incoming,
                1,
            )),
            builder.add_vertex(DiagramVertex::external_in_connection(
                "out0",
                2,
                ExternalState::Outgoing,
                0,
            )),
            builder.add_vertex(DiagramVertex::external_in_connection(
                "out1",
                3,
                ExternalState::Outgoing,
                1,
            )),
        ];
        let interaction = model.vertex_rule_id("V").unwrap();
        let left = builder.add_vertex(DiagramVertex::interaction("left", interaction));
        let right = builder.add_vertex(DiagramVertex::interaction("right", interaction));
        let marker = builder.add_vertex(DiagramVertex::interaction("marker", interaction));
        let targets = match attachment {
            LeftRightAttachment::Identity => [left, left, right, right],
            LeftRightAttachment::GlobalSwap => [right, right, left, left],
            LeftRightAttachment::PartialSwap => [right, left, left, right],
        };
        for (index, (&external, &target)) in external.iter().zip(&targets).enumerate() {
            let (source, target) = if index < 2 {
                (external, target)
            } else {
                (target, external)
            };
            builder
                .add_edge(source, target, DiagramEdge::new(particle(model, 25), false))
                .unwrap();
        }
        builder
            .add_edge(left, right, DiagramEdge::new(particle(model, 25), false))
            .unwrap();
        builder
            .add_edge(left, marker, DiagramEdge::new(particle(model, 25), false))
            .unwrap();
        builder.build().unwrap()
    }

    fn exact_options() -> GraphGroupingOptions {
        GraphGroupingOptions {
            check_canonical_numerator: true,
            ..GraphGroupingOptions::default()
        }
    }

    #[test]
    fn groups_identical_numerators() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x+y", 25, 25),
                diagram(&model, "g1", "y+x", 25, 25),
            ],
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            false,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].master, 0);
        assert_eq!(grouped.groups[0].members[1].ratio, test_atom("1"));
        assert!(
            grouped.diagrams[0]
                .overall_factor()
                .get_all_symbols(true)
                .contains(&symbol!(
                    "feynkit_generator_factor::NumeratorDependentGrouping"
                ))
        );
    }

    #[test]
    fn lmb_routing_and_projector_local_ids_do_not_split_numerator_groups() {
        let model = Arc::new(model());
        let external_first = routed_projector_diagram(&model, "g0", false);
        let internal_first = routed_projector_diagram(&model, "g1", true);
        assert_ne!(
            external_first.loop_momentum_basis().loop_edges,
            internal_first.loop_momentum_basis().loop_edges
        );
        assert_ne!(external_first.projector(), internal_first.projector());
        let grouped = group_diagrams(
            vec![external_first, internal_first],
            &model,
            &NumeratorGrouping::Identical(GraphGroupingOptions {
                check_canonical_numerator: false,
                number_of_numerical_samples: 2,
                ..GraphGroupingOptions::default()
            }),
            false,
        )
        .unwrap();

        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members.len(), 2);
    }

    #[test]
    fn groups_numerators_up_to_sign() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x+y", 25, 25),
                diagram(&model, "g1", "-x-y", 25, 25),
            ],
            &model,
            &NumeratorGrouping::UpToSign(exact_options()),
            false,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members[1].ratio, test_atom("-1"));
    }

    #[test]
    fn scalar_ratio_is_member_over_master() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x+y", 25, 25),
                diagram(&model, "g1", "2*x+2*y", 25, 25),
            ],
            &model,
            &NumeratorGrouping::UpToScalar(exact_options()),
            false,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members[1].ratio, test_atom("2"));
    }

    #[test]
    fn exact_scalar_ratio_preserves_symbolic_color_dimension() {
        let model = Arc::new(model());
        let nc = Atom::var(CS.nc);
        let common_tensor = symbol!("FeynKit::grouping_test_tensor").call(1);
        let master = &nc * &common_tensor;
        let candidate = (nc.clone().pow(Atom::num(2)) - Atom::one()) * &common_tensor;
        let ratio = compare_atoms(
            &candidate,
            &master,
            ComparisonMode::UpToScalar,
            &scalar_names(&model),
        )
        .unwrap();

        assert_eq!(
            ratio,
            ((nc.clone().pow(Atom::num(2)) - Atom::one()) / nc).cancel()
        );

        let other_tensor = symbol!("FeynKit::other_grouping_test_tensor").call(1);
        assert!(
            compare_atoms(
                &common_tensor,
                &other_tensor,
                ComparisonMode::UpToScalar,
                &scalar_names(&model),
            )
            .is_none(),
            "a quotient retaining tensor heads must not be accepted as a scalar ratio"
        );
    }

    #[test]
    fn topology_uses_internal_mass_and_spin_but_exact_external_species() {
        let model = Arc::new(model());
        let diagrams = vec![
            diagram(&model, "g0", "1", 25, 25),
            diagram(&model, "g1", "1", 35, 25),
            diagram(&model, "g2", "1", 45, 25),
            diagram(&model, "g3", "1", 25, 35),
        ];
        let grouped = group_diagrams(
            diagrams.clone(),
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            false,
        )
        .unwrap();
        assert_eq!(
            grouped
                .groups
                .iter()
                .map(|group| {
                    group
                        .members
                        .iter()
                        .map(|member| member.source_diagram)
                        .collect()
                })
                .collect::<Vec<Vec<_>>>(),
            vec![vec![0, 1], vec![2], vec![3]]
        );

        let exact_species = GraphGroupingOptions {
            differentiate_particle_masses_only: false,
            ..exact_options()
        };
        let grouped = group_diagrams(
            diagrams,
            &model,
            &NumeratorGrouping::Identical(exact_species),
            false,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 4);
    }

    #[test]
    fn distinct_cut_metadata_is_merged_into_the_denominator_master() {
        let model = Arc::new(model());
        let diagrams = vec![cut_diagram(&model, false), cut_diagram(&model, true)];
        let source_cuts = diagrams
            .iter()
            .flat_map(|diagram| diagram.cuts().iter().cloned())
            .collect::<BTreeSet<_>>();
        let grouped = group_diagrams(
            diagrams,
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            false,
        )
        .unwrap();

        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members.len(), 2);
        assert_eq!(grouped.diagrams.len(), 1);
        grouped.diagrams[0].validate().unwrap();
        assert_eq!(
            grouped.diagrams[0]
                .cuts()
                .iter()
                .cloned()
                .collect::<BTreeSet<_>>(),
            source_cuts,
            "every source cut must survive on the grouped denominator master"
        );
    }

    #[test]
    fn grouped_partitions_are_remapped_from_permuted_member_ids() {
        let model = Arc::new(model());
        let master = partitioned_diagram(&model, "master", false, false);
        let member = partitioned_diagram(&model, "member", true, false);
        assert_ne!(
            master.cuts()[0].left.half_edges,
            member.cuts()[0].left.half_edges,
            "the regression requires different stable half-edge frames"
        );
        let expected_cuts = master.cuts().to_vec();
        let expected_candidates = master.topology_threshold_candidates().to_vec();

        let grouped = group_diagrams(
            vec![master, member],
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            false,
        )
        .unwrap();

        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members.len(), 2);
        assert_eq!(grouped.diagrams[0].cuts(), expected_cuts);
        assert_eq!(
            grouped.diagrams[0].topology_threshold_candidates(),
            expected_candidates
        );
        grouped.diagrams[0].validate().unwrap();
    }

    #[test]
    fn grouped_partitions_preserve_source_side_across_global_left_right_swap() {
        let model = Arc::new(model());
        let master = partitioned_diagram(&model, "master", false, false);
        let member = partitioned_diagram(&model, "member", true, true);
        let expected_cuts = master.cuts().to_vec();
        let expected_candidates = master.topology_threshold_candidates().to_vec();

        let grouped = group_diagrams(
            vec![master, member],
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            true,
        )
        .unwrap();

        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members.len(), 2);
        assert_eq!(grouped.diagrams[0].cuts(), expected_cuts);
        assert_eq!(
            grouped.diagrams[0].topology_threshold_candidates(),
            expected_candidates,
            "the globally swapped member must not add an opposite-shift candidate"
        );
        grouped.diagrams[0].validate().unwrap();
    }

    #[test]
    fn only_one_global_left_right_swap_is_in_the_grouping_orbit() {
        let model = Arc::new(model());
        let identity = left_right_diagram(&model, LeftRightAttachment::Identity);
        let global = left_right_diagram(&model, LeftRightAttachment::GlobalSwap);
        let partial = left_right_diagram(&model, LeftRightAttachment::PartialSwap);

        let grouped = group_diagrams(
            vec![identity.clone(), global],
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            true,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members.len(), 2);

        let grouped = group_diagrams(
            vec![identity, partial],
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            true,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 2);
    }

    #[test]
    fn removes_only_proven_zero_numerators_and_preserves_source_ids() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x-x", 25, 25),
                diagram(&model, "g1", "1", 25, 25),
            ],
            &model,
            &NumeratorGrouping::OnlyDetectZeroes,
            false,
        )
        .unwrap();
        assert_eq!(grouped.zero_numerator_count, 1);
        assert_eq!(grouped.diagrams.len(), 1);
        assert_eq!(grouped.groups[0].members[0].diagram, 0);
        assert_eq!(grouped.groups[0].members[0].source_diagram, 1);
        assert_eq!(grouped.diagrams[0].name(), "g1");
    }

    #[test]
    fn zero_projectors_are_removed_before_grouping() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![diagram_with_projector(&model, "zero", "1", "0", 25, 25)],
            &model,
            &NumeratorGrouping::OnlyDetectZeroes,
            false,
        )
        .unwrap();

        assert_eq!(grouped.zero_numerator_count, 1);
        assert!(grouped.diagrams.is_empty());
        assert!(grouped.groups.is_empty());
    }

    #[test]
    fn coupling_definitions_are_applied_before_zero_detection() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![diagram(&model, "coupling-zero", "UFO::GC-1", 25, 25)],
            &model,
            &NumeratorGrouping::OnlyDetectZeroes,
            false,
        )
        .unwrap();

        assert_eq!(grouped.zero_numerator_count, 1);
        assert!(grouped.diagrams.is_empty());
    }
}
