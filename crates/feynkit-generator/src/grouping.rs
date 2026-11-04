use crate::generation::FACTOR_HEADS;
#[cfg(test)]
use crate::generation::cp_merge_diagnostic;
use feynkit_kinematics::symbols as kinematics_symbols;
#[cfg(test)]
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use feynkit_graph::expressions::evaluate_overall_factor;
use feynkit_graph::symbols;
use feynkit_graph::{
    DiagramEndpoint, DiagramError, DiagramHalfEdge, EdgeId, ExternalLeg, ExternalState,
    FeynmanDiagram, VertexId,
};
use feynkit_model::{Model, ModelError, ParameterId, ParticleId};
use idenso::color::{CS, ColorSimplifySettings};
use idenso::dirac::{AGS, spinor_matrix_structure};
use idenso::rep_symbols::RS;
use idenso::representations::{Bispinor, ColorAdjoint, ColorFundamental, ColorSextet};
use idenso::tensor::inference::{InterfaceInference, TensorInferenceError};
use idenso::tensor::{AlgebraContraction, AlgebraSettings, ContractSettings, SymbolicTensor};
use idenso::{CookMode, CookSettings, cof, color_idx};
use linnet::half_edge::builder::HedgeGraphBuilder;
use linnet::half_edge::canon::{CanonicalForm, ignore};
use linnet::half_edge::involution::{Flow, Orientation};
use linnet::half_edge::{HedgeGraph, NoData, NodeIndex};
use spenso::network::library::LibraryTensor;
use spenso::network::library::function_lib::{INBUILTS, Panic};
use spenso::network::library::symbolic::{ExplicitKey, TensorLibrary};
use spenso::network::parsing::{
    ParseSettings as TensorParseSettings, ShadowedStructure, ShorthandParsing,
};
use spenso::network::store::NetworkStore;
use spenso::network::{MinResultRank, Network, Sequential};
use spenso::shadowing::{Collectable, TensorCollectFilter};
use spenso::structure::Canonicalized;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
use spenso::structure::slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind, SlotError};
use spenso::tensors::data::{DataTensor, StorageTensor};
use spenso::tensors::parametric::ParamTensor;
use spenso_hep_lib::{
    CHARGE_CONJUGATION_WEYL_COMPONENTS, gamma_data_weyl, gamma5_weyl_data, proj_m_data_weyl,
    proj_p_data_weyl,
};
use symbolica::atom::representation::FunView;
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol};
use symbolica::id::{AliasedAtom, Replacement};
#[cfg(test)]
use symbolica::parser::ParseSettings;
use symbolica::{function, symbol};
use thiserror::Error;

use crate::generation::DiagramPair;
use crate::{DiagramGroup, GenerationError, GraphGroupingOptions, GroupMember, Grouping, Run};

#[derive(Debug, Error)]
pub enum GroupingError {
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error("failed to construct a topology grouping key: {0}")]
    TopologyConstruction(String),
    #[error("source diagram index {0} does not fit in a symbolic integer")]
    SourceDiagramIndexOverflow(usize),
    #[error(
        "failed to merge physical cuts from group members {members:?} into master diagram {master:?}: {source}"
    )]
    CutMerge {
        master: Option<String>,
        members: Vec<Option<String>>,
        #[source]
        source: Box<DiagramError>,
    },
    #[error(
        "failed to merge topology threshold candidates from group members {members:?} into master diagram {master:?}: {source}"
    )]
    ThresholdCandidateMerge {
        master: Option<String>,
        members: Vec<Option<String>>,
        #[source]
        source: Box<DiagramError>,
    },
    #[error(
        "failed to transport grouped partition from member diagram {member:?} to master diagram {master:?}: {message}"
    )]
    PartitionTransport {
        member: Option<String>,
        master: Option<String>,
        message: String,
    },
    #[error("failed to compare numerators of diagrams {candidate:?} and {master:?}: {source}")]
    NumeratorComparison {
        candidate: Option<String>,
        master: Option<String>,
        #[source]
        source: Box<TensorInferenceError>,
    },
    #[error("failed to simplify numerator color for diagram {diagram:?}: {message}")]
    ColorSimplification {
        diagram: Option<String>,
        message: String,
    },
    #[error(
        "failed to evaluate numerator tensors for diagram {diagram:?}, sample {sample}: {message}"
    )]
    TensorEvaluation {
        diagram: Option<String>,
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
    Hedge(usize, usize),
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
            Self::Hedge(owner, local) => write!(formatter, "hedge({owner},{local})"),
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
    #[error("{name} has an unsupported number of non-negative integer arguments: {actual}")]
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
    if !(1..=2).contains(&function.get_nargs()) {
        return Err(GroupingIndexError::InvalidArity {
            name,
            actual: function.get_nargs(),
        });
    }
    let mut arguments = function.iter();
    let owner = natural_index(arguments.next().expect("arity checked"))?;
    let local = arguments
        .next()
        .map(natural_index)
        .transpose()?
        .unwrap_or(0);
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
                if head == feynkit_graph::symbols::hedge_index() {
                    owned_index(function, "gammalooprs::hedge", Self::Hedge)
                } else if head == feynkit_graph::symbols::vertex_index() {
                    owned_index(function, "gammalooprs::vertex", Self::Vertex)
                } else if head == feynkit_graph::symbols::edge_index() {
                    owned_index(function, "gammalooprs::edge", Self::Edge)
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
            GroupingIndex::Hedge(owner, local) => {
                AbstractIndex::Named(symbols::hedge_index().into(), owner, local).to_atom()
            }
            GroupingIndex::Vertex(owner, local) => {
                AbstractIndex::Named(symbols::vertex_index().into(), owner, local).to_atom()
            }
            GroupingIndex::Edge(owner, local) => {
                AbstractIndex::Named(symbols::edge_index().into(), owner, local).to_atom()
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

impl GroupingOutcome {
    /// Singleton groups that move each physical owner to the output.
    fn ungrouped(
        diagrams: Vec<DiagramPair>,
        source_diagrams: &[usize],
        zero_numerator_count: usize,
    ) -> Result<Self, GenerationError> {
        let mut masters = Vec::with_capacity(diagrams.len());
        let mut groups = Vec::with_capacity(diagrams.len());
        for (index, (pair, &source_diagram)) in
            diagrams.into_iter().zip(source_diagrams).enumerate()
        {
            groups.push(DiagramGroup::singleton(
                index,
                source_diagram,
                None,
                pair.comparison().overall_factor().clone(),
            ));
            masters.push(pair.into_representative()?);
        }
        Ok(Self {
            diagrams: masters,
            groups,
            zero_numerator_count,
        })
    }
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
        external: ExternalLeg,
    },
}

type TopologyKey = CanonicalForm<NoData, TopologyEdge, NoData>;

struct PreparedNumerator {
    exact: SymbolicTensor<PartialStructure>,
    sample_source: SymbolicTensor<PartialStructure>,
    // Lorentz components have been evaluated, but color ports can remain open.
    samples: Vec<SymbolicTensor<PartialStructure>>,
}

impl PreparedNumerator {
    // This bounds a complete interpolation certificate, not a statistical
    // zero test. Unsupported syntax or a larger grid stays unproved.
    const MAX_COEFFICIENT_PROBES: usize = 256;

    fn sample_degrees(value: AtomView<'_>) -> Option<BTreeMap<Symbol, usize>> {
        let (product, parts) = match value {
            AtomView::Num(number) => {
                return matches!(
                    number.get_coeff_view(),
                    symbolica::coefficient::CoefficientView::Natural(..)
                        | symbolica::coefficient::CoefficientView::Large(..)
                )
                .then(BTreeMap::new);
            }
            AtomView::Var(variable) => return Some(BTreeMap::from([(variable.get_symbol(), 1)])),
            AtomView::Pow(power) => {
                let exponent = usize::try_from(power.get_exp()).ok()?;
                let mut degrees = Self::sample_degrees(power.get_base())?;
                for degree in degrees.values_mut() {
                    *degree = degree.checked_mul(exponent)?;
                }
                return degrees
                    .values()
                    .try_fold(1usize, |size, degree| {
                        size.checked_mul(degree.checked_add(1)?)
                    })
                    .filter(|size| *size <= Self::MAX_COEFFICIENT_PROBES)
                    .map(|_| degrees);
            }
            AtomView::Add(sum) => (false, sum.iter().collect::<Vec<_>>()),
            AtomView::Mul(factors) => (true, factors.iter().collect::<Vec<_>>()),
            _ => return None,
        };
        let mut degrees = BTreeMap::<Symbol, usize>::new();
        for part in parts {
            for (symbol, degree) in Self::sample_degrees(part)? {
                let current = degrees.entry(symbol).or_default();
                *current = if product {
                    current.checked_add(degree)?
                } else {
                    (*current).max(degree)
                };
            }
            let size = degrees.values().try_fold(1usize, |size, degree| {
                size.checked_mul(degree.checked_add(1)?)
            })?;
            if size > Self::MAX_COEFFICIENT_PROBES {
                return None;
            }
        }
        Some(degrees)
    }

    fn finite_samples_equal(left: &Atom, right: &Atom) -> bool {
        let difference = (left - right).collect_factors();
        if difference.is_zero() {
            return true;
        }
        let factors = match difference.as_view() {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            other => vec![other],
        };
        factors.into_iter().any(|factor| {
            let Some(degrees) = Self::sample_degrees(factor) else {
                return false;
            };
            if degrees.is_empty() {
                return factor.is_zero();
            }
            let count = degrees.values().map(|degree| degree + 1).product::<usize>();
            // Degree d in each variable requires every point of the Cartesian
            // (d+1)-grid. Exact Gaussian-rational arithmetic certifies zero;
            // substituting coordinates never distributes the factored body.
            (0..count).all(|mut point| {
                let replacements = degrees
                    .iter()
                    .map(|(symbol, degree)| {
                        let coordinate = point % (degree + 1);
                        point /= degree + 1;
                        Replacement::new(Atom::var(*symbol), Atom::num(coordinate))
                    })
                    .collect::<Vec<_>>();
                factor.replace_multiple(&replacements).is_zero()
            })
        })
    }
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

/// What the caller knows about each input diagram of [`group_diagrams`]: its interference
/// class (groups never mix classes) and its index in the caller's diagram list. Without them,
/// all diagrams share one class and keep their input positions.
#[derive(Clone, Copy, Default)]
pub(crate) struct DiagramLabels<'a> {
    pub(crate) classes: Option<&'a [usize]>,
    pub(crate) source_indices: Option<&'a [usize]>,
}

pub(crate) fn group_diagrams(
    diagrams: Vec<DiagramPair>,
    model: &Model,
    grouping: &Grouping,
    symmetrize_left_right: bool,
    generation_options: &Run,
    filter_zero_color: bool,
    labels: DiagramLabels<'_>,
) -> Result<GroupingOutcome, GenerationError> {
    let DiagramLabels {
        classes,
        source_indices,
    } = labels;
    let source_index = |index: usize| source_indices.map_or(index, |indices| indices[index]);
    if generation_options.cancellation_requested()
        || (matches!(grouping, Grouping::None) && !filter_zero_color)
    {
        let source_diagrams = (0..diagrams.len()).map(source_index).collect::<Vec<_>>();
        return GroupingOutcome::ungrouped(diagrams, &source_diagrams, 0);
    }

    let total = diagrams.len();
    // Numerator comparisons reuse the same color reduction; enabling the
    // standalone filter never adds a second color-algebra pass.
    let preparation_stage = if filter_zero_color {
        crate::Stage::ZeroColor
    } else {
        crate::Stage::GroupingPreparation
    };
    generation_options.report_progress(preparation_stage, 0usize, Some(total));
    let compares_numerators = matches!(
        grouping,
        Grouping::Identical(_) | Grouping::UpToSign(_) | Grouping::UpToScalar(_)
    );
    let (compares_canonical_numerators, concrete_canonical_color) = match grouping {
        Grouping::Identical(options)
        | Grouping::UpToSign(options)
        | Grouping::UpToScalar(options) => (
            options.check_canonical_numerator,
            options.fully_numerical_substitution,
        ),
        Grouping::None => (false, false),
    };
    let symmetric_polarizations = match grouping {
        Grouping::Identical(options)
        | Grouping::UpToSign(options)
        | Grouping::UpToScalar(options) => options.symmetric_polarizations,
        Grouping::None => false,
    };
    let progress = Mutex::new((0usize, 0usize));
    let prepare = |(source_diagram, pair): (usize, DiagramPair)| {
        if generation_options.cancellation_requested() {
            // Keep unchecked inputs on cancellation. Their preparation is
            // unused: the incomplete result below contains singleton groups.
            return Ok(Some((
                pair,
                source_diagram,
                PreparedNumerator {
                    exact: SymbolicTensor::infer(Atom::one()).expect("literal scalar one"),
                    sample_source: SymbolicTensor::infer(Atom::one()).expect("literal scalar one"),
                    samples: Vec::new(),
                },
            )));
        }
        let diagram = pair.comparison();
        #[cfg(test)]
        cp_merge_diagnostic::record(|| {
            json!({
                "stage": "comparison_source", "source": source_index(source_diagram),
                "diagram": cp_merge_diagnostic::operator(diagram),
            })
        });
        let color_error = |error: TensorInferenceError| GroupingError::ColorSimplification {
            diagram: diagram.name().map(str::to_owned),
            message: error.to_string(),
        };
        let complete = diagram
            .numerator()
            .map_err(|error| GroupingError::ColorSimplification {
                diagram: diagram.name().map(str::to_owned),
                message: error.to_string(),
            })?
            .multiply(diagram.numerator_prefactor())
            .map_err(|error| color_error(error.into()))?
            .multiply(diagram.projector())
            .map_err(|error| color_error(error.into()))?;
        let complete = complete
            .with_checked_expression(model.expand_couplings(complete.expression()))
            .map_err(&color_error)?;
        let color_reps = [
            ColorAdjoint {}.into(),
            ColorFundamental {}.into(),
            ColorSextet {}.into(),
        ];
        let complete = if compares_numerators {
            complete
                .map_representation_dimensions(|rep, dimension| {
                    if rep.base() == (ColorAdjoint {}).into() {
                        Atom::var(CS.nc).pow(2) - 1
                    } else if rep.base() == (ColorFundamental {}).into() {
                        Atom::var(CS.nc)
                    } else {
                        dimension.to_owned()
                    }
                })
                .map_err(&color_error)?
        } else {
            complete
        };
        let simplified = complete
            .simplify_algebra(&AlgebraSettings {
                color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                contract: if compares_canonical_numerators {
                    AlgebraContraction::Fully
                } else if compares_numerators {
                    // Numerical Lorentz contraction belongs after momentum sampling.
                    AlgebraContraction::Selected(color_reps.to_vec())
                } else {
                    AlgebraContraction::None
                },
                ..Default::default()
            })
            .map_err(&color_error)?;
        let is_zero = if compares_numerators {
            simplified
                .coefficients_are_zero(TensorCollectFilter::Reps(color_reps))
                .map_err(&color_error)?
                .is_true()
        } else {
            // Retain the existing conservative exact-zero policy. Color-only
            // algebra keeps unselected tensor factors in their native owner.
            simplified
                .expression()
                .collect_factors()
                .check_zero_with(
                    &symbolica::atom::ZeroTestSettings::default()
                        .finite_field_samples(0)
                        .numerical_samples(0)
                        .exact_fallback(false),
                )
                .is_true()
        };
        #[cfg(test)]
        cp_merge_diagnostic::record(|| {
            json!({
                "stage": "numerator_admission", "source": source_index(source_diagram), "zero": is_zero,
            })
        });
        if is_zero {
            return Ok(None);
        }
        let sample_source = if compares_numerators {
            normalize_momentum_routing(diagram, &simplified).map_err(&color_error)?
        } else {
            SymbolicTensor::infer(Atom::one()).expect("literal scalar one")
        };
        // Canonical comparison and numerical sampling retain separate typed
        // values, preserving the existing optional canonicalization policy.
        let exact = if compares_canonical_numerators {
            let exact_source = if symmetric_polarizations {
                normalize_symmetric_polarizations(&sample_source).map_err(&color_error)?
            } else {
                sample_source.clone()
            };
            let exact_source = if concrete_canonical_color {
                let concrete = exact_source
                    .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))
                    .map_err(&color_error)?;
                concrete
                    .with_checked_expression(
                        concrete
                            .expression()
                            .replace(CS.nc)
                            .with(3)
                            .replace(color_idx!(2, cof!(RS.d_)))
                            .with(Atom::num(1) / Atom::num(2)),
                    )
                    .map_err(&color_error)?
            } else {
                exact_source
            };
            // Factoring can change which tensor sum anchors the dummy frame.
            // Finish scalar normalization before assigning canonical indices.
            let exact_source = exact_source
                .with_algebra_result(exact_source.expression().replace_map_bottom_up(
                    |value, _, out| {
                        if matches!(value, AtomView::Add(_)) {
                            **out = value.collect_numeric_content();
                        }
                    },
                ))
                .map_err(&color_error)?;
            let exact_source = if matches!(grouping, Grouping::UpToScalar(_)) {
                exact_source
                    .with_algebra_result(exact_source.expression().collect_factors())
                    .map_err(|error| tensor_evaluation_error(diagram, 0, error))?
            } else {
                exact_source
            };
            exact_source.canonize().map_err(|error| {
                tensor_evaluation_error(diagram, 0, format!("exact canonicalization: {error}"))
            })?
        } else {
            SymbolicTensor::infer(Atom::one()).expect("literal scalar one")
        };
        #[cfg(test)]
        cp_merge_diagnostic::record(|| {
            json!({
                "stage": "canonical_numerator", "source": source_index(source_diagram),
                "exact": exact.expression().to_canonical_string(),
            })
        });
        Ok(Some((
            pair,
            source_diagram,
            PreparedNumerator {
                exact,
                sample_source,
                samples: Vec::new(),
            },
        )))
    };
    let preparations = generation_options
        .map_ordered(diagrams.into_iter().enumerate().collect(), |input| {
            let result = prepare(input);
            // Serialize progress callbacks while preserving input order in the
            // indexed collection, independently of worker completion order.
            let mut progress = progress.lock().unwrap();
            progress.0 += 1;
            progress.1 += usize::from(matches!(result, Ok(None)));
            generation_options.report_progress_with_zero_count(
                preparation_stage,
                progress.0,
                Some(total),
                progress.1,
            );
            result
        })
        .into_iter()
        .collect::<Result<Vec<_>, GroupingError>>()?;
    let mut retained = Vec::with_capacity(total);
    let mut source_diagrams = Vec::with_capacity(total);
    let mut prepared = Vec::with_capacity(total);
    for (diagram, source_diagram, numerator) in preparations.into_iter().flatten() {
        retained.push(diagram);
        source_diagrams.push(source_diagram);
        prepared.push(numerator);
    }
    let mut zero_numerator_count = total - retained.len();

    let grouping = if generation_options.cancellation_requested() {
        &Grouping::None
    } else {
        grouping
    };
    let (options, mode) = match grouping {
        Grouping::Identical(options) => (options, ComparisonMode::Identical),
        Grouping::UpToSign(options) => (options, ComparisonMode::UpToSign),
        Grouping::UpToScalar(options) => (options, ComparisonMode::UpToScalar),
        Grouping::None => {
            return GroupingOutcome::ungrouped(
                retained,
                &source_diagrams
                    .iter()
                    .copied()
                    .map(source_index)
                    .collect::<Vec<_>>(),
                zero_numerator_count,
            );
        }
    };
    let scalar_names = scalar_names(model);
    let total = retained.len();
    generation_options.report_progress_with_zero_count(
        crate::Stage::GroupingSamples,
        0usize,
        Some(total),
        zero_numerator_count,
    );
    let mut buckets = BTreeMap::<(Option<usize>, TopologyKey), Vec<usize>>::new();
    let mut canonical_frames = Vec::with_capacity(retained.len());
    for (index, (pair, numerator)) in retained.iter().zip(&mut prepared).enumerate() {
        let diagram = pair.comparison();
        let canonical = canonical_topology(diagram, model, options, symmetrize_left_right)?;
        let key = canonical.key.clone();
        canonical_frames.push(canonical);
        #[cfg(test)]
        cp_merge_diagnostic::record(|| {
            json!({
                "stage": "group_frame", "source": source_index(source_diagrams[index]),
                "vertices": canonical_frames[index].vertex_map,
                "swapped": canonical_frames[index].left_right_swapped,
            })
        });
        numerator.samples = numerical_tensor_samples(&numerator.sample_source, diagram, options)?;
        let mut zero_samples = !numerator.samples.is_empty();
        for (sample, value) in numerator.samples.iter().enumerate() {
            if !PreparedNumerator::finite_samples_equal(value.expression(), &Atom::Zero)
                && !value
                    .coefficients_are_zero(TensorCollectFilter::Reps([
                        ColorAdjoint {}.into(),
                        ColorFundamental {}.into(),
                        ColorSextet {}.into(),
                    ]))
                    .map(|proof| proof.is_true())
                    .map_err(|error| tensor_evaluation_error(diagram, sample, error))?
            {
                zero_samples = false;
                break;
            }
        }
        if zero_samples {
            zero_numerator_count += 1;
        } else {
            buckets
                .entry((classes.map(|classes| classes[source_diagrams[index]]), key))
                .or_default()
                .push(index);
        }
        generation_options.report_progress_with_zero_count(
            crate::Stage::GroupingSamples,
            index + 1,
            Some(total),
            zero_numerator_count,
        );
    }
    let total = buckets.values().map(Vec::len).sum();
    let mut completed = 0usize;
    generation_options.report_progress_with_zero_count(
        crate::Stage::GroupingComparison,
        completed,
        Some(total),
        zero_numerator_count,
    );
    let mut groups = Vec::new();
    for indices in buckets.into_values() {
        let mut bucket_groups: Vec<DiagramGroup> = Vec::new();
        for diagram in indices {
            let mut matching_group = None;
            for (group_index, group) in bucket_groups.iter().enumerate() {
                if let Some(ratio) = compare(
                    &prepared[diagram],
                    &prepared[group.master],
                    mode,
                    options,
                    &scalar_names,
                )
                .map_err(|source| GroupingError::NumeratorComparison {
                    candidate: retained[diagram].comparison().name().map(str::to_owned),
                    master: retained[group.master]
                        .comparison()
                        .name()
                        .map(str::to_owned),
                    source: Box::new(source),
                })? {
                    matching_group = Some((group_index, ratio));
                    break;
                }
            }
            if let Some((group_index, ratio)) = matching_group {
                bucket_groups[group_index].members.push(GroupMember {
                    source_diagram: source_index(source_diagrams[diagram]),
                    source_identity: None,
                    diagram,
                    ratio,
                    overall_factor: retained[diagram].comparison().overall_factor().clone(),
                });
            } else {
                bucket_groups.push(DiagramGroup {
                    master: diagram,
                    members: vec![GroupMember {
                        source_diagram: source_index(source_diagrams[diagram]),
                        source_identity: None,
                        diagram,
                        ratio: Atom::one(),
                        overall_factor: retained[diagram].comparison().overall_factor().clone(),
                    }],
                });
            }
            completed += 1;
            generation_options.report_progress_with_zero_count(
                crate::Stage::GroupingComparison,
                completed,
                Some(total),
                zero_numerator_count,
            );
        }
        groups.extend(bucket_groups);
    }
    // Compare complete scalar weights after numerator ratios have been proved.
    // Opposite members can cancel even though neither numerator is zero. Keep
    // the graph numerators and retained provenance of surviving groups intact.
    groups.retain(|group| {
        let weight = Atom::add_many(group.members.iter().map(|member| {
            &member.ratio * evaluate_overall_factor(member.overall_factor.as_view())
        }));
        // These are scalar normalization weights, not graph numerators.
        // Exact rescaling can leave distinct denominators or algebraic
        // constants such as sqrt(1/2)*sqrt(2). No sampled zero is accepted.
        let keep = !weight
            .check_zero_with(
                &symbolica::atom::ZeroTestSettings::default()
                    .finite_field_samples(0)
                    .numerical_samples(0)
                    .exact_fallback(true),
            )
            .is_true();
        #[cfg(test)]
        cp_merge_diagnostic::record(|| {
            json!({
                "stage": "group_component", "keep": keep,
                "master_source": source_index(source_diagrams[group.master]),
                "weight": weight.to_canonical_string(),
                "members": group.members.iter().map(|member| json!({
                    "source": member.source_diagram, "ratio": member.ratio.to_canonical_string(),
                    "factor": member.overall_factor.to_canonical_string(),
                })).collect::<Vec<_>>(),
            })
        });
        keep
    });
    groups.sort_by_key(|group| group.master);

    generation_options.report_progress_with_zero_count(
        crate::Stage::Grouping,
        0usize,
        None,
        zero_numerator_count,
    );
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
fn normalize_momentum_routing(
    diagram: &FeynmanDiagram,
    source: &SymbolicTensor<PartialStructure>,
) -> Result<SymbolicTensor<PartialStructure>, TensorInferenceError> {
    let mut numerator = source.expression().clone();
    let basis = diagram.loop_momentum_basis();
    let external_connections = diagram
        .edges()
        .filter_map(|(id, _, edge)| edge.external.as_ref().map(|leg| (id, leg.connection)))
        .collect::<BTreeMap<_, _>>();
    numerator = basis.route_expression(&numerator);
    let args = symbol!("feynkit_generator::routing_args___");
    let replacements = basis
        .external_edges
        .iter()
        .enumerate()
        .filter_map(|(slot, edge)| {
            external_connections.get(edge).map(|connection| {
                Replacement::new(
                    function!(kinematics_symbols::external_momentum(), slot, args).to_pattern(),
                    function!(kinematics_symbols::external_momentum(), *connection, args)
                        .to_pattern(),
                )
            })
        })
        .collect::<Vec<_>>();
    numerator = numerator.replace_multiple(&replacements);
    let polarization_symbols = [
        feynkit_graph::symbols::u(),
        feynkit_graph::symbols::ubar(),
        feynkit_graph::symbols::v(),
        feynkit_graph::symbols::vbar(),
        feynkit_graph::symbols::epsilon(),
        feynkit_graph::symbols::epsilonbar(),
    ];

    numerator = numerator.replace_map(|term, _, output| {
        let AtomView::Fun(function) = term else {
            return;
        };
        if polarization_symbols.contains(&function.get_symbol()) && function.get_nargs() == 2 {
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
    source.with_checked_expression(numerator)
}

/// Give left/right partners the same symbolic polarization tensor for exact
/// numerator comparisons, matching the shared numerical samples selected by
/// `symmetric_polarizations`.
fn normalize_symmetric_polarizations(
    numerator: &SymbolicTensor<PartialStructure>,
) -> Result<SymbolicTensor<PartialStructure>, TensorInferenceError> {
    numerator.with_checked_expression(numerator.expression().replace_map(|term, _, output| {
        let AtomView::Fun(function) = term else {
            return;
        };
        let normalized_head = if [
            feynkit_graph::symbols::epsilon(),
            feynkit_graph::symbols::epsilonbar(),
        ]
        .contains(&function.get_symbol())
        {
            feynkit_graph::symbols::epsilon()
        } else if [
            feynkit_graph::symbols::u(),
            feynkit_graph::symbols::ubar(),
            feynkit_graph::symbols::v(),
            feynkit_graph::symbols::vbar(),
        ]
        .contains(&function.get_symbol())
        {
            feynkit_graph::symbols::u()
        } else {
            return;
        };
        **output = FunctionBuilder::new(normalized_head)
            .add_args(function.iter())
            .finish();
    }))
}

fn momentum_edge(argument: AtomView<'_>) -> Option<EdgeId> {
    i64::try_from(argument)
        .ok()
        .and_then(|edge| usize::try_from(edge).ok())
        .map(EdgeId)
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
                    member: diagram.name().map(str::to_owned),
                    master: master.name().map(str::to_owned),
                    message: format!("unknown source edge {}", half_edge.edge.0),
                }
            })?;
            match half_edge.endpoint {
                DiagramEndpoint::Source => endpoints.source,
                DiagramEndpoint::Target => endpoints.target,
            }
            .ok_or_else(|| GroupingError::PartitionTransport {
                member: diagram.name().map(str::to_owned),
                master: master.name().map(str::to_owned),
                message: format!("half-edge {:?} is absent", half_edge),
            })
        })
        .collect()
}

/// Grouping compares denominator topologies, but half-edge IDs belong to each input graph.
/// Transport the member's vertex partition through the canonical topology frame before
/// rebuilding its half-edges on the master.
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
            member: member.name().map(str::to_owned),
            master: master.name().map(str::to_owned),
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
                member: member.name().map(str::to_owned),
                master: master.name().map(str::to_owned),
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
                    member: member.name().map(str::to_owned),
                    master: master.name().map(str::to_owned),
                    message: format!("source vertex {} has no canonical image", vertex.0),
                }
            })?;
            canonical_to_master.get(canonical).copied().ok_or_else(|| {
                GroupingError::PartitionTransport {
                    member: member.name().map(str::to_owned),
                    master: master.name().map(str::to_owned),
                    message: format!("canonical vertex {canonical} has no master preimage"),
                }
            })
        })
        .collect::<Result<BTreeSet<_>, _>>()?;

    Ok(master
        .edges()
        .flat_map(|(edge, endpoints, _)| {
            [
                endpoints
                    .source
                    .is_some_and(|vertex| mapped_vertices.contains(&vertex))
                    .then_some(DiagramHalfEdge {
                        edge,
                        endpoint: DiagramEndpoint::Source,
                    }),
                endpoints
                    .target
                    .is_some_and(|vertex| mapped_vertices.contains(&vertex))
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
    diagrams: Vec<DiagramPair>,
    canonical_frames: &[CanonicalTopology],
    groups: Vec<DiagramGroup>,
) -> Result<(Vec<FeynmanDiagram>, Vec<DiagramGroup>), GenerationError> {
    // Compute merged partitions while all inputs are borrowed. Afterwards each sorted master
    // has one consumer; nonmasters can be dropped.
    let plans = groups
        .into_iter()
        .map(|group| {
            let topology_master = diagrams[group.master].comparison();
            let master_frame = &canonical_frames[group.master];
            let mut cut_partitions = Vec::new();
            let mut threshold_partitions = Vec::new();
            let mut members = Vec::with_capacity(group.members.len());
            for member in &group.members {
                let source = diagrams[member.diagram].comparison();
                let source_frame = &canonical_frames[member.diagram];
                members.push(source.name().map(str::to_owned));
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
            cut_partitions.sort();
            cut_partitions.dedup();
            threshold_partitions.sort();
            threshold_partitions.dedup();
            let factor = if group.members.len() == 1 {
                topology_master.overall_factor().clone()
            } else {
                let mut factor = Atom::Zero;
                for member in &group.members {
                    let source_diagram = i64::try_from(member.source_diagram).map_err(|_| {
                        GroupingError::SourceDiagramIndexOverflow(member.source_diagram)
                    })?;
                    factor += function!(
                        FACTOR_HEADS.numerator_dependent_grouping,
                        Atom::num(source_diagram),
                        member.ratio.clone(),
                        member.overall_factor.clone()
                    );
                }
                factor
            };
            Ok::<_, GenerationError>((group, cut_partitions, threshold_partitions, factor, members))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut masters = Vec::with_capacity(plans.len());
    let mut collapsed_groups = Vec::with_capacity(plans.len());
    let mut inputs = diagrams.into_iter().enumerate();
    for (group, cut_partitions, threshold_partitions, factor, members) in plans {
        let (_, pair) = inputs
            .find(|(index, _)| *index == group.master)
            .expect("group masters are unique and sorted by input index");
        let merge = |mut master: FeynmanDiagram| -> Result<_, GroupingError> {
            if !cut_partitions.is_empty() || !master.cuts().is_empty() {
                let master_name = master.name().map(str::to_owned);
                master = master
                    .with_cut_partitions(cut_partitions)
                    .map_err(|source| GroupingError::CutMerge {
                        master: master_name,
                        members: members.clone(),
                        source: Box::new(source),
                    })?;
            }
            if !threshold_partitions.is_empty()
                || !master.topology_threshold_candidates().is_empty()
            {
                let master_name = master.name().map(str::to_owned);
                master = master
                    .with_topology_threshold_partitions(threshold_partitions)
                    .map_err(|source| GroupingError::ThresholdCandidateMerge {
                        master: master_name,
                        members,
                        source: Box::new(source),
                    })?;
            }
            Ok(master.with_overall_factor(factor))
        };
        let pair = match pair {
            DiagramPair::Shared(master) => DiagramPair::Shared(merge(master)?),
            DiagramPair::Distinct {
                comparison,
                representative,
                reversed_edges,
            } => DiagramPair::Distinct {
                comparison: merge(comparison)?,
                representative,
                reversed_edges,
            },
        };
        let output_index = masters.len();
        masters.push(pair.into_representative()?);
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
    numerator: &SymbolicTensor<PartialStructure>,
    diagram: &FeynmanDiagram,
    options: &GraphGroupingOptions,
) -> Result<Vec<SymbolicTensor<PartialStructure>>, GroupingError> {
    (0..options.number_of_numerical_samples)
        .map(|sample| evaluate_tensor_sample(numerator, diagram, options, sample))
        .collect()
}

fn evaluate_tensor_sample(
    numerator: &SymbolicTensor<PartialStructure>,
    diagram: &FeynmanDiagram,
    options: &GraphGroupingOptions,
    sample: usize,
) -> Result<SymbolicTensor<PartialStructure>, GroupingError> {
    if numerator.is_zero() {
        let color_reps = [
            ColorAdjoint {}.into(),
            ColorFundamental {}.into(),
            ColorSextet {}.into(),
        ];
        if numerator
            .structure()
            .logical_slots()
            .iter()
            .any(|slot| !color_reps.contains(&slot.rep().rep.base()))
        {
            return Err(tensor_evaluation_error(
                diagram,
                sample,
                "Lorentz sampling retained open tensor ports",
            ));
        }
        return numerator
            .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))
            .map_err(|error| tensor_evaluation_error(diagram, sample, error));
    }
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
        Replacement::new(color_idx!(2, cof!(RS.d_)), Atom::num(1) / Atom::num(2)),
    ];
    let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
    let sectors = numerator
        .coefficient_list(TensorCollectFilter::Reps([
            ColorAdjoint {}.into(),
            ColorFundamental {}.into(),
            ColorSextet {}.into(),
        ]))
        .map_err(|error| {
            tensor_evaluation_error(diagram, sample, format!("color-sector separation: {error}"))
        })?;
    let contributions = sectors
        .into_iter()
        .map(|(color, lorentz)| {
            // Separating a color coefficient can expose a closed metric loop
            // that the mixed tensor expression retained. Complete that isolated
            // contraction before comparing it with a colorless scalar sample.
            let color = color
                .contract(ContractSettings::default())
                .map_err(|error| {
                    tensor_evaluation_error(diagram, sample, format!("color contraction: {error}"))
                })?;
            let lowered: Canonicalized<GroupingTensorNetwork> = lorentz
                .into_network(&library, &functions, &settings, Some(&cooking))
                .map_err(|error| {
                    tensor_evaluation_error(diagram, sample, format!("Lorentz parsing: {error}"))
                })?;
            let (mut network, layout) = lowered.into_parts();
            if layout.order() != 0 {
                return Err(tensor_evaluation_error(
                    diagram,
                    sample,
                    "Lorentz sampling retained open tensor ports",
                ));
            }
            for scalar in &mut network.store.scalar {
                *scalar = scalar.replace_multiple(&scalar_replacements);
            }
            network
                .execute::<Sequential, MinResultRank, _, _, _>(&library, &functions)
                .map_err(|error| {
                    tensor_evaluation_error(diagram, sample, format!("Lorentz execution: {error}"))
                })?;
            let mut scalar: Atom = network
                .result_scalar()
                .map_err(|error| {
                    tensor_evaluation_error(diagram, sample, format!("scalar extraction: {error}"))
                })?
                .into();
            // Contraction leaves exact finite component coefficients in different
            // nested Horner forms. Collect those coefficients before comparing
            // samples; retain every function and composite power as an opaque
            // factor, so this cannot distribute a graph numerator or denominator.
            let mut exact = true;
            scalar.visitor(&mut |part| {
                if let AtomView::Num(number) = part {
                    exact &= matches!(
                        number.get_coeff_view(),
                        symbolica::coefficient::CoefficientView::Natural(..)
                            | symbolica::coefficient::CoefficientView::Large(..)
                    );
                }
                exact
            });
            if exact {
                let mut used = scalar.get_all_symbols(true);
                let mut aliases = BTreeMap::<Atom, Atom>::new();
                let mut serial = 0usize;
                let mut collected = scalar.replace_map(|part, _, output| {
                    let opaque = match part {
                        AtomView::Fun(_) => true,
                        AtomView::Pow(power) => {
                            !matches!(power.get_base(), AtomView::Var(_))
                                || !i64::try_from(power.get_exp())
                                    .is_ok_and(|exponent| exponent > 0)
                        }
                        _ => false,
                    };
                    if opaque {
                        let alias = aliases.entry(part.to_owned()).or_insert_with(|| {
                            loop {
                                let symbol = symbol!(&format!(
                                    "feynkit_generator::sample_coefficient_{serial}"
                                ));
                                serial += 1;
                                if used.insert(symbol) {
                                    break Atom::var(symbol);
                                }
                            }
                        });
                        **output = alias.clone();
                    }
                });
                let normalize_content = |value: &Atom| {
                    value.replace_map_bottom_up(|part, _, output| {
                        if matches!(part, AtomView::Add(_)) {
                            **output = part.collect_numeric_content();
                        }
                    })
                };
                collected = normalize_content(&collected);
                loop {
                    let next =
                        normalize_content(&collected.expand_num().collect_horner::<Symbol>(None))
                            .collect_factors();
                    if next.as_view().get_byte_size() >= collected.as_view().get_byte_size() {
                        break;
                    }
                    collected = next;
                }
                let mut collected = AliasedAtom::from(collected);
                for (original, alias) in aliases {
                    collected.register_alias(alias, original);
                }
                scalar = collected.into_inner();
            }
            let color = color.canonize().map_err(|error| {
                tensor_evaluation_error(diagram, sample, format!("color canonicalization: {error}"))
            })?;
            let color = color
                .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))
                .map_err(|error| tensor_evaluation_error(diagram, sample, error))?;
            color
                .with_checked_expression(
                    (color.expression() * scalar).replace_multiple(&color_replacements),
                )
                .map_err(|error| tensor_evaluation_error(diagram, sample, error))
        })
        .collect::<Result<Vec<_>, GroupingError>>()?;
    let mut contributions = contributions.into_iter();
    let Some(mut sum) = contributions.next() else {
        return numerator
            .with_checked_expression(Atom::Zero)
            .map_err(|error| tensor_evaluation_error(diagram, sample, error));
    };
    for contribution in contributions {
        sum = sum
            .try_add(&contribution)
            .map_err(|error| tensor_evaluation_error(diagram, sample, error))?;
    }
    Ok(sum)
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
            kinematics_symbols::loop_momentum(),
            Atom::num(loop_index),
            &format!("loop-momentum:{loop_index}"),
            false,
            options,
            sample,
        )?;
    }

    let mut external_legs = diagram
        .edges()
        .filter_map(|(_, endpoints, edge)| {
            edge.external
                .as_ref()
                .map(|leg| (leg.clone(), edge.particle, endpoints))
        })
        .collect::<Vec<_>>();
    external_legs.sort_by_key(|(external, particle, _)| {
        (external.connection, external.state, particle.index())
    });

    let mut inserted_momenta = BTreeSet::new();
    for (external, particle_id, endpoints) in external_legs {
        if inserted_momenta.insert(external.connection) {
            insert_sample_vector(
                &mut library,
                kinematics_symbols::external_momentum(),
                Atom::num(external.connection),
                &format!("external-momentum:{}", external.connection),
                false,
                options,
                sample,
            )?;
        }

        let particle = diagram.model().particle_by_id(particle_id)?;
        // A sewn carrier owns both endpoint wavefunctions, just as in the
        // finalized projector; an unsewn leg retains its declared state.
        for (vertex, state) in [
            (endpoints.source, ExternalState::Outgoing),
            (endpoints.target, ExternalState::Incoming),
        ] {
            if vertex.is_none() {
                continue;
            }
            let state = if endpoints.source.is_some() && endpoints.target.is_some() {
                state
            } else {
                external.state
            };
            let head = match (particle.spin, state, particle.is_antiparticle()) {
                (2, ExternalState::Incoming, false) => feynkit_graph::symbols::u(),
                (2, ExternalState::Incoming, true) => feynkit_graph::symbols::vbar(),
                (2, ExternalState::Outgoing, false) => feynkit_graph::symbols::ubar(),
                (2, ExternalState::Outgoing, true) => feynkit_graph::symbols::v(),
                (3, ExternalState::Incoming, _) => feynkit_graph::symbols::epsilon(),
                (3, ExternalState::Outgoing, _) => feynkit_graph::symbols::epsilonbar(),
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
    }
    Ok(library)
}

fn insert_dirac_tensors(library: &mut GroupingTensorLibrary) {
    let gamma = gamma_data_weyl(
        AGS.gamma_strct::<GroupingIndex>(4),
        Atom::num(1),
        Atom::num(0),
    )
    .map_canonical(|tensor| {
        ParamTensor::composite(DataTensor::Sparse(
            tensor.map_data(|value| value.re + Atom::i() * value.im),
        ))
    });
    library.insert_explicit(gamma);
    let gamma5 = gamma5_weyl_data(
        AGS.gamma5_strct::<GroupingIndex>(4),
        Atom::num(1),
        Atom::num(0),
    )
    .map_canonical(|tensor| {
        ParamTensor::composite(DataTensor::Sparse(
            tensor.map_data(|value| value.re + Atom::i() * value.im),
        ))
    });
    library.insert_explicit(gamma5);
    let projm = proj_m_data_weyl(
        AGS.projm_strct::<GroupingIndex>(4),
        Atom::num(1),
        Atom::num(0),
    )
    .map_canonical(|tensor| {
        ParamTensor::composite(DataTensor::Sparse(
            tensor.map_data(|value| value.re + Atom::i() * value.im),
        ))
    });
    library.insert_explicit(projm);
    let projp = proj_p_data_weyl(
        AGS.projp_strct::<GroupingIndex>(4),
        Atom::num(1),
        Atom::num(0),
    )
    .map_canonical(|tensor| {
        ParamTensor::composite(DataTensor::Sparse(
            tensor.map_data(|value| value.re + Atom::i() * value.im),
        ))
    });
    library.insert_explicit(projp);
    library
        .insert_explicit_sparse(
            spinor_matrix_structure::<GroupingIndex>(AGS.charge_conjugation, 4),
            CHARGE_CONJUGATION_WEYL_COMPONENTS
                .map(|(indices, sign)| (indices.to_vec(), Atom::num(sign))),
            Atom::Zero,
        )
        .expect("charge-conjugation entries have valid four-dimensional spinor indices");
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
    let tensor = ParamTensor::from_dense(key.into_canonical(), data).map_err(|error| {
        GroupingError::TensorEvaluation {
            diagram: Some("sample library".to_owned()),
            sample,
            message: error.to_string(),
        }
    })?;
    library.insert_explicit(Canonicalized::identity(tensor));
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
    let tensor = ParamTensor::from_dense(key.into_canonical(), data).map_err(|error| {
        GroupingError::TensorEvaluation {
            diagram: Some("sample library".to_owned()),
            sample,
            message: error.to_string(),
        }
    })?;
    library.insert_explicit(Canonicalized::identity(tensor));
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
        diagram: diagram.name().map(str::to_owned),
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
) -> Result<Option<Atom>, TensorInferenceError> {
    // All provided completed samples belong to the certificate, including when
    // canonical equality would otherwise return before inspecting them.
    if candidate
        .samples
        .iter()
        .chain(&master.samples)
        .any(|value| !value.expression().is_finite())
    {
        return Ok(None);
    }
    if options.check_canonical_numerator
        && let Some(ratio) = compare_atoms(&candidate.exact, &master.exact, mode, scalar_names)?
    {
        // A contracted representation dimension can survive as a scalar
        // coefficient. Its decoded Nc becomes visible only at this quotient
        // boundary, so apply the requested concrete-color choice here too.
        let ratio = if options.fully_numerical_substitution {
            ratio.replace(CS.nc).with(3)
        } else {
            ratio
        };
        return Ok(ratio.is_finite().then_some(ratio));
    }
    let mut first: Option<Atom> = None;
    for (candidate, master) in candidate.samples.iter().zip(&master.samples) {
        if !InterfaceInference::additive_interfaces_match(candidate.structure(), master.structure())
        {
            return Ok(None);
        }
        // This certificate is restricted to completed finite component samples.
        // The optional raw tensor comparison above keeps its existing policy.
        let ratio =
            if PreparedNumerator::finite_samples_equal(candidate.expression(), master.expression())
            {
                Some(Atom::one())
            } else if matches!(mode, ComparisonMode::UpToSign)
                && PreparedNumerator::finite_samples_equal(
                    candidate.expression(),
                    &(-master.expression()),
                )
            {
                Some(-Atom::one())
            } else {
                compare_atoms(candidate, master, mode, scalar_names)?
            };
        let Some(ratio) = ratio else {
            return Ok(None);
        };
        if let Some(first) = &first {
            // These ratios contain only model/color scalars: all tensor heads
            // have already been excluded. Compare their rational normalizations
            // without distributing either source numerator.
            if ratio != *first
                && !(&ratio - first).together().is_zero()
                && !PreparedNumerator::finite_samples_equal(
                    candidate.expression(),
                    &(first * master.expression()),
                )
            {
                return Ok(None);
            }
        } else {
            first = Some(ratio);
        }
    }
    Ok(first)
}

fn compare_atoms(
    candidate: &SymbolicTensor<PartialStructure>,
    master: &SymbolicTensor<PartialStructure>,
    mode: ComparisonMode,
    scalar_names: &BTreeSet<String>,
) -> Result<Option<Atom>, TensorInferenceError> {
    if !InterfaceInference::additive_interfaces_match(candidate.structure(), master.structure()) {
        return Ok(None);
    }
    // Reject both operands before equality or division can conceal undefined
    // values (equal infinities or a finite value divided by infinity).
    if !candidate.expression().is_finite() || !master.expression().is_finite() {
        return Ok(None);
    }
    if expressions_equal(candidate, master)? {
        return Ok(Some(Atom::num(1)));
    }
    if matches!(mode, ComparisonMode::UpToSign) && expressions_equal(candidate, &(-master.clone()))?
    {
        return Ok(Some(Atom::num(-1)));
    }
    if !matches!(mode, ComparisonMode::UpToScalar) || master.is_zero() {
        return Ok(None);
    }
    let ratio = (candidate.expression() / master.expression()).collect_factors();
    // Closed color contractions can emit an encoded representation dimension
    // as a scalar weight. Decode only that metadata in the cancelled quotient;
    // keep both tensor bodies and encoded index identities untouched.
    let ratio = ratio
        .replace_map(|view, _, out| {
            if let AtomView::Var(variable) = view
                && let Some(dimension) = CookSettings::dimension_payload(variable.get_symbol())
            {
                **out = dimension;
            }
        })
        .collect_factors();
    if expression_is_scalar(&ratio, scalar_names) {
        return Ok(ratio.is_finite().then_some(ratio));
    }
    // Completed component samples can place exact numerical content differently
    // in equivalent colour-sector sums. Pull content only from existing sums
    // in the quotient, as canonical preparation already does. This neither
    // distributes a numerator nor changes either admitted operand or its history.
    let ratio = ratio
        .replace_map_bottom_up(|part, _, out| {
            if matches!(part, AtomView::Add(_)) {
                **out = part.collect_numeric_content();
            }
        })
        .collect_factors();
    // Polynomial conversion of a graph numerator is deliberately excluded.
    // Dummy indices can be alpha-equivalent across separately prepared graphs;
    // accept an exact quotient only after every tensor head has disappeared.
    Ok((expression_is_scalar(&ratio, scalar_names) && ratio.is_finite()).then_some(ratio))
}

fn expressions_equal(
    left: &SymbolicTensor<PartialStructure>,
    right: &SymbolicTensor<PartialStructure>,
) -> Result<bool, TensorInferenceError> {
    let filter = TensorCollectFilter::Reps([
        ColorAdjoint {}.into(),
        ColorFundamental {}.into(),
        ColorSextet {}.into(),
    ]);
    if right.is_zero() {
        return left
            .coefficients_are_zero(filter)
            .map(|proof| proof.is_true());
    }
    if left.expression().is_zero() {
        return right
            .coefficients_are_zero(filter)
            .map(|proof| proof.is_true());
    }
    left.coefficients_equal(right, filter)
        .map(|proof| proof.is_true())
}

fn expression_is_scalar(expression: &Atom, scalar_names: &BTreeSet<String>) -> bool {
    let mut scalar = true;
    expression.visitor(&mut |atom| {
        if !scalar {
            return false;
        }
        let symbol = match atom {
            AtomView::Var(variable) => variable.get_symbol(),
            AtomView::Fun(function)
                if [Symbol::CONJ, INBUILTS.conj].contains(&function.get_symbol()) =>
            {
                // Conjugation broadcasts over both scalars and tensors. Its
                // argument must still pass this check; the head alone does
                // not certify a scalar ratio (notably for CKM couplings).
                scalar = function.get_nargs() == 1;
                return scalar;
            }
            AtomView::Fun(function) => function.get_symbol(),
            _ => return true,
        };
        // Scalar functions own their metadata: representation labels inside
        // a Casimir or Dynkin index are not open tensor factors.
        if symbol.is_scalar() {
            return false;
        }
        scalar = symbol == CS.nc || scalar_names.contains(symbol.get_stripped_name());
        scalar
    });
    scalar
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
    let mut topologies = vec![topology_key_variant(diagram, model, options, false)?];
    if symmetrize_left_right && is_sewn(diagram) {
        topologies.push(topology_key_variant(diagram, model, options, true)?);
    }
    Ok(topologies)
}

/// Whether every external edge of `diagram` is a sewn initial-state edge.
fn is_sewn(diagram: &FeynmanDiagram) -> bool {
    let mut external = diagram
        .edges()
        .filter(|(_, _, edge)| edge.external.is_some())
        .peekable();
    external.peek().is_some()
        && external
            .all(|(_, endpoints, _)| endpoints.source.is_some() && endpoints.target.is_some())
}

/// The canonical topology of `diagram`, with the sewn initial-state edges reversed for the
/// global left/right exchange when `left_right_swapped`.
fn topology_key_variant(
    diagram: &FeynmanDiagram,
    model: &Model,
    options: &GraphGroupingOptions,
    left_right_swapped: bool,
) -> Result<CanonicalTopology, GroupingError> {
    let mut graph = HedgeGraphBuilder::new();
    let nodes: Vec<NodeIndex> = diagram.vertices().map(|_| graph.add_node(())).collect();
    for (_, endpoints, edge) in diagram.edges() {
        let (orientation, color) = if let Some(external) = &edge.external {
            // External edges retain both their exact particle species and
            // direction: exchanging different mass-degenerate external states
            // would change the process rather than merely its topology.
            (
                Orientation::Default,
                TopologyEdge::External {
                    particle: edge.particle,
                    directed: edge.directed,
                    external: external.clone(),
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
            (Orientation::Undirected, TopologyEdge::Internal(color))
        };
        let (mut source, mut target) = (endpoints.source, endpoints.target);
        if left_right_swapped && edge.external.is_some() {
            std::mem::swap(&mut source, &mut target);
        }
        match (source, target) {
            (Some(source), Some(target)) => {
                graph.add_edge(nodes[source.0], nodes[target.0], color, orientation)
            }
            (Some(_), None) | (None, Some(_)) if edge.external.is_none() => {
                return Err(GroupingError::TopologyConstruction(
                    "dangling edge has no external metadata".into(),
                ));
            }
            (Some(vertex), None) => {
                graph.add_external_edge(nodes[vertex.0], color, orientation, Flow::Source)
            }
            (None, Some(vertex)) => {
                graph.add_external_edge(nodes[vertex.0], color, orientation, Flow::Sink)
            }
            (None, None) => {
                return Err(GroupingError::TopologyConstruction(
                    "edge has no attached vertex".into(),
                ));
            }
        }
    }
    let graph: HedgeGraph<TopologyEdge, ()> = graph.build();
    let canonical = graph.canonize_by(ignore, Clone::clone, ignore);
    Ok(CanonicalTopology {
        vertex_map: nodes
            .into_iter()
            .map(|node| canonical.canonical_node(node).0)
            .collect(),
        key: canonical.into_canonical_form(),
        left_right_swapped,
    })
}

#[cfg(test)]
mod tests {
    use crate::FeynmanDiagramRulesExt;
    use crate::Process;
    use crate::{GenerationControl, Stage};
    use feynkit_graph as graph_api;
    use feynkit_graph::NumeratorInterface;
    use feynkit_model::ParameterCard;
    use idenso::CookSourceFilter;
    use idenso::IndexTooling;
    use serde::Deserialize;
    use serde_json::Value;
    use spenso::iterators::IteratableTensor;
    use spenso::network::library::symbolic::ETS;
    use spenso::structure::TensorDataLayout;
    use spenso::structure::TensorStructure;
    use spenso_hep_lib as hep_lib;
    use std::env;
    use std::fs;
    use std::path::PathBuf;
    use std::process::Command;
    use std::sync::Arc;
    use std::sync::atomic::AtomicBool;
    use std::thread;
    use symbolica::atom::ZeroTestSettings;

    use feynkit_graph::{
        DiagramCut, DiagramCutSide, DiagramEdge, DiagramEndpoint, DiagramError, DiagramHalfEdge,
        DiagramVertex, ExternalState, VertexId, VertexSlot,
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

    // C gamma^T C^-1 = -gamma. In a closed spinor circuit the inserted
    // C matrices cancel; the spinor metric is unchanged. Keep each matrix
    // product and propagator sum factorized, as in the sewn self-energy proof.
    fn conjugate_spinor_circuits(input: &Atom) -> Atom {
        input.replace_map(|atom, _, output| {
            if let AtomView::Fun(function) = atom
                && function.get_symbol() == AGS.gamma
            {
                assert_eq!(function.get_nargs(), 3);
                **output = -FunctionBuilder::new(AGS.gamma)
                    .add_args([function.get(1), function.get(0), function.get(2)])
                    .finish();
            }
        })
    }

    #[test]
    #[ignore = "bounded proof over frozen scalar observations; requires FEYNKIT_ORACLE_DUMP"]
    fn scalar_changed_factorized_components() {
        let directory = PathBuf::from(env::var("FEYNKIT_ORACLE_DUMP").unwrap());
        let observations = ["baseline", "candidate"].map(|label| {
            let data = fs::read(directory.join(format!("{label}.json"))).unwrap();
            serde_json::from_slice::<Value>(&data).unwrap()
        });
        let data = fs::read(directory.join("components.json")).unwrap();
        let components: Value = serde_json::from_slice(&data).unwrap();
        let model = Model::qcd_qed();
        let options = GraphGroupingOptions {
            number_of_numerical_samples: 0,
            fully_numerical_substitution: true,
            check_canonical_numerator: true,
            symmetric_polarizations: true,
            ..GraphGroupingOptions::default()
        };
        let scalars = scalar_names(&model);
        let mut results = Vec::new();
        for component in components.as_array().unwrap() {
            let anchor = component["candidate_groups"][0]["master_source"]
                .as_u64()
                .unwrap();
            let exact = |side: usize, source: u64| {
                let event = observations[side]["events"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .find(|event| {
                        event["stage"] == "canonical_numerator" && event["source"] == source
                    })
                    .unwrap();
                let expression = test_atom(event["exact"].as_str().unwrap());
                SymbolicTensor::infer(expression).unwrap()
            };
            let master = exact(1, anchor);
            let master = PreparedNumerator {
                sample_source: master.clone(),
                exact: master,
                samples: Vec::new(),
            };
            let mut totals = Vec::new();
            let mut certificates = Vec::new();
            let mut complete = true;
            for (side, label) in ["baseline", "candidate"].into_iter().enumerate() {
                let mut total = Atom::Zero;
                for group in component[format!("{label}_groups")].as_array().unwrap() {
                    for member in group["members"].as_array().unwrap() {
                        let source = member["source"].as_u64().unwrap();
                        let numerator = exact(side, source);
                        let prepared = PreparedNumerator {
                            sample_source: numerator.clone(),
                            exact: numerator,
                            samples: Vec::new(),
                        };
                        let mut ratio = compare(
                            &prepared,
                            &master,
                            ComparisonMode::UpToScalar,
                            &options,
                            &scalars,
                        )
                        .unwrap();
                        let mut charge_conjugated = false;
                        if ratio.is_none() {
                            assert_eq!(
                                prepared
                                    .exact
                                    .structure()
                                    .canonical()
                                    .external_structure_iter()
                                    .count(),
                                0,
                                "charge conjugation requires closed spinor circuits"
                            );
                            let conjugated = conjugate_spinor_circuits(prepared.exact.expression());
                            assert_eq!(
                                conjugate_spinor_circuits(&conjugated),
                                *prepared.exact.expression()
                            );
                            let exact = prepared
                                .exact
                                .with_checked_expression(conjugated)
                                .unwrap()
                                .canonize()
                                .unwrap();
                            let conjugated = PreparedNumerator {
                                sample_source: exact.clone(),
                                exact,
                                samples: Vec::new(),
                            };
                            ratio = compare(
                                &conjugated,
                                &master,
                                ComparisonMode::UpToScalar,
                                &options,
                                &scalars,
                            )
                            .unwrap();
                            charge_conjugated = true;
                        }
                        if let Some(ratio) = &ratio {
                            let factor = test_atom(member["factor"].as_str().unwrap());
                            total += ratio * evaluate_overall_factor(factor.as_view());
                        } else {
                            complete = false;
                        }
                        certificates.push(json!({
                            "side": label, "source": source,
                            "charge_conjugated_closed_spinor_circuits": charge_conjugated,
                            "ratio_to_anchor": ratio.map(|ratio| ratio.to_canonical_string()),
                        }));
                    }
                }
                totals.push(total);
            }
            let difference = (&totals[0] - &totals[1]).collect_factors();
            let equal = complete
                && difference
                    .check_zero_with(
                        &ZeroTestSettings::default()
                            .finite_field_samples(0)
                            .numerical_samples(0)
                            .exact_fallback(true),
                    )
                    .is_true();
            results.push(json!({
                "component": component["component"], "anchor_candidate_source": anchor,
                "certificates": certificates, "complete": complete, "equal": equal,
                "difference": difference.to_canonical_string(),
            }));
        }
        let output = directory.join("proof-result.json");
        assert!(!output.exists(), "proof output must be new");
        fs::write(&output, serde_json::to_vec(&results).unwrap()).unwrap();
        assert!(results.iter().all(|result| result["equal"] == true));
    }

    // Authored fixture expressions enter at the source boundary. These exact
    // oracles intentionally compare decoded expressions after the production
    // path has retained its carrier through all intermediate operations.
    fn evaluate_tensor_sample(
        source: &Atom,
        diagram: &FeynmanDiagram,
        options: &GraphGroupingOptions,
        sample: usize,
    ) -> Result<Atom, GroupingError> {
        let source = feynkit_graph::NumeratorInterface::admit(source.clone(), None)
            .map_err(|error| tensor_evaluation_error(diagram, sample, error))?;
        let result = super::evaluate_tensor_sample(&source, diagram, options, sample)?;
        Ok(feynkit_graph::NumeratorInterface::export(&result).0)
    }

    fn expressions_equal(left: &Atom, right: &Atom) -> Result<bool, TensorInferenceError> {
        let admit = |source: &Atom| {
            feynkit_graph::NumeratorInterface::admit(source.clone(), None)
                .map_err(|error| TensorInferenceError::Invalid(error.to_string()))
        };
        super::expressions_equal(&admit(left)?, &admit(right)?)
    }

    fn compare_atoms(
        left: &Atom,
        right: &Atom,
        mode: ComparisonMode,
        scalars: &BTreeSet<String>,
    ) -> Result<Option<Atom>, TensorInferenceError> {
        let admit = |source: &Atom| {
            feynkit_graph::NumeratorInterface::admit(source.clone(), None)
                .map_err(|error| TensorInferenceError::Invalid(error.to_string()))
        };
        super::compare_atoms(&admit(left)?, &admit(right)?, mode, scalars)
    }

    fn normalize_momentum_routing(diagram: &FeynmanDiagram, source: Atom) -> Atom {
        let source = feynkit_graph::NumeratorInterface::admit(source, None).unwrap();
        let result = super::normalize_momentum_routing(diagram, &source).unwrap();
        feynkit_graph::NumeratorInterface::export(&result).0
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

    fn leg(name: &str, index: usize, state: ExternalState, connection: usize) -> ExternalLeg {
        ExternalLeg {
            name: name.to_owned(),
            index,
            state,
            connection,
        }
    }

    /// Builds fixtures with the endpoint slots in edge insertion order, as `add_edge` assigns
    /// them, for internal edges and generated legs alike.
    struct FixtureBuilder {
        builder: feynkit_graph::FeynmanDiagramBuilder,
        degrees: BTreeMap<VertexId, usize>,
    }

    impl FixtureBuilder {
        fn new(builder: feynkit_graph::FeynmanDiagramBuilder) -> Self {
            Self {
                builder,
                degrees: BTreeMap::new(),
            }
        }

        fn slot(&mut self, vertex: VertexId) -> VertexSlot {
            let degree = self.degrees.entry(vertex).or_default();
            *degree += 1;
            VertexSlot(*degree - 1)
        }

        fn add_vertex(&mut self, vertex: DiagramVertex) -> VertexId {
            self.builder.add_vertex(vertex)
        }

        fn add_edge(&mut self, source: VertexId, target: VertexId, edge: DiagramEdge) -> EdgeId {
            let (source_slot, target_slot) = (self.slot(source), self.slot(target));
            self.builder
                .add_edge_with_slots(source, target, edge, source_slot, target_slot)
                .unwrap()
        }

        fn add_leg(
            &mut self,
            vertex: VertexId,
            flow: Flow,
            edge: DiagramEdge,
            leg: ExternalLeg,
        ) -> EdgeId {
            let slot = self.slot(vertex);
            self.builder
                .add_generation_external(vertex, flow, edge, slot, leg)
                .unwrap()
        }
    }

    fn diagram_with_atoms(
        model: &Arc<Model>,
        name: &str,
        numerator: Atom,
        projector: Atom,
        internal_pdg: i64,
        external_pdg: i64,
    ) -> FeynmanDiagram {
        let mut builder = FixtureBuilder::new(
            FeynmanDiagram::builder(Arc::clone(model), name)
                .numerator(feynkit_graph::NumeratorInterface::admit(numerator, None).unwrap())
                .projector(feynkit_graph::NumeratorInterface::admit(projector, None).unwrap()),
        );
        let interaction = model.vertex_rule_id("V").unwrap();
        let left = builder.add_vertex(DiagramVertex::interaction("left", interaction));
        let right = builder.add_vertex(DiagramVertex::interaction("right", interaction));
        builder.add_leg(
            left,
            Flow::Sink,
            DiagramEdge::new(particle(model, external_pdg), false),
            leg("in", 0, ExternalState::Incoming, 0),
        );
        builder.add_edge(
            left,
            right,
            DiagramEdge::new(particle(model, internal_pdg), false),
        );
        builder.add_leg(
            right,
            Flow::Source,
            DiagramEdge::new(particle(model, external_pdg), false),
            leg("out", 1, ExternalState::Outgoing, 1),
        );
        builder.builder.build().unwrap()
    }

    fn routed_projector_diagram(
        model: &Arc<Model>,
        name: &str,
        internal_edges_first: bool,
    ) -> FeynmanDiagram {
        let mut builder = FixtureBuilder::new(FeynmanDiagram::builder(Arc::clone(model), name));
        let interaction = model.vertex_rule_id("V").unwrap();
        let left = builder.add_vertex(DiagramVertex::interaction("left", interaction));
        let right = builder.add_vertex(DiagramVertex::interaction("right", interaction));
        let scalar = || DiagramEdge::new(particle(model, 25), false);
        let vector = || DiagramEdge::new(particle(model, 21), false);
        let add_internals = |builder: &mut FixtureBuilder| {
            let first = builder.add_edge(left, right, scalar());
            let second = builder.add_edge(left, right, scalar());
            if internal_edges_first { first } else { second }
        };
        let add_incoming = |builder: &mut FixtureBuilder| {
            builder.add_leg(
                left,
                Flow::Sink,
                vector(),
                leg("in", 0, ExternalState::Incoming, 0),
            )
        };
        let (incoming_edge, loop_edge) = if internal_edges_first {
            let loop_edge = add_internals(&mut builder);
            (add_incoming(&mut builder), loop_edge)
        } else {
            let incoming_edge = add_incoming(&mut builder);
            (incoming_edge, add_internals(&mut builder))
        };
        builder.add_leg(
            right,
            Flow::Source,
            vector(),
            leg("out", 1, ExternalState::Outgoing, 1),
        );
        let builder = builder.builder;
        let endpoint = feynkit_graph::symbols::hedge_index().call((
            (2 * incoming_edge.0 + 1) as i64,
            if internal_edges_first { 2_i64 } else { 1_i64 },
        ));
        let index = Minkowski {}.new_rep(4).to_symbolic([endpoint]);
        let numerator = FunctionBuilder::new(kinematics_symbols::momentum())
            .add_arg(loop_edge.0 as i64)
            .add_arg(index.clone())
            .finish();
        let projector = FunctionBuilder::new(feynkit_graph::symbols::epsilon())
            .add_arg(incoming_edge.0 as i64)
            .add_arg(index)
            .finish();
        let diagram = builder
            .numerator(feynkit_graph::NumeratorInterface::admit(numerator, None).unwrap())
            .projector(feynkit_graph::NumeratorInterface::admit(projector, None).unwrap())
            .build()
            .unwrap();
        let AtomView::Fun(momentum) = diagram.numerator().unwrap().expression().as_view() else {
            panic!("the fixture numerator is one edge momentum");
        };
        let loop_edge = EdgeId(usize::try_from(momentum.iter().next().unwrap()).unwrap());
        let diagram = diagram.with_loop_momentum_edges(&[loop_edge]).unwrap();
        assert_eq!(diagram.loop_momentum_basis().loop_edges, vec![loop_edge]);
        diagram
    }

    fn cut_diagram(model: &Arc<Model>, mirrored: bool) -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), "cut-line");
        let vertices = (0..3)
            .map(|index| {
                builder.add_vertex(DiagramVertex {
                    name: format!("interaction{index}"),
                    interaction: None,
                    numerator: SymbolicTensor::infer(Atom::one()).unwrap(),
                })
            })
            .collect::<Vec<_>>();
        let mut external = DiagramEdge::new(particle(model, 25), false);
        external.external = Some(ExternalLeg {
            name: "incoming".into(),
            index: 0,
            state: ExternalState::Incoming,
            connection: 0,
        });
        builder
            .add_edge(vertices[0], vertices[2], external)
            .unwrap();
        for pair in vertices.windows(2) {
            builder
                .add_edge(
                    pair[0],
                    pair[1],
                    DiagramEdge::new(particle(model, 25), false),
                )
                .unwrap();
        }
        let diagram = builder.build().unwrap();
        // The global transformation exchanges both physical sides at once;
        // the two distinct cuts cross the first or second internal propagator.
        let selected = if mirrored {
            &vertices[1..]
        } else {
            &vertices[2..]
        };
        let (left, right): (Vec<_>, Vec<_>) = diagram.half_edges().partition(|half_edge| {
            let node = diagram
                .underlying()
                .node_id(diagram.half_edge_id(*half_edge).unwrap());
            selected.contains(&VertexId(node.0))
        });
        let diagram = diagram.with_cut_partitions(vec![(left, right)]).unwrap();
        diagram.validate().unwrap();
        diagram
    }

    fn partitioned_diagram(
        model: &Arc<Model>,
        name: &str,
        permuted_ids: bool,
        global_left_right_swap: bool,
    ) -> FeynmanDiagram {
        /// An edge end: a vertex or the far end of a generated leg.
        #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
        enum End {
            Vertex(VertexId),
            Incoming,
            Outgoing,
        }
        let mut builder = FixtureBuilder::new(FeynmanDiagram::builder(Arc::clone(model), name));
        let internal = |name: &str| DiagramVertex {
            name: name.to_owned(),
            interaction: None,

            numerator: SymbolicTensor::infer(Atom::one()).unwrap(),
        };
        let (left, right, marker) = if permuted_ids {
            let marker = builder.add_vertex(internal("marker"));
            let right = builder.add_vertex(internal("right"));
            let left = builder.add_vertex(internal("left"));
            (left, right, marker)
        } else {
            let left = builder.add_vertex(internal("left"));
            let right = builder.add_vertex(internal("right"));
            let marker = builder.add_vertex(internal("marker"));
            (left, right, marker)
        };
        let (source_side, target_side) = if global_left_right_swap {
            (right, left)
        } else {
            (left, right)
        };
        let scalar = || DiagramEdge::new(particle(model, 25), false);
        let mut edges = Vec::new();
        let mut add_edge = |source, target| {
            let edge = match (source, target) {
                (End::Vertex(source), End::Vertex(target)) => {
                    builder.add_edge(source, target, scalar())
                }
                (End::Incoming, End::Vertex(target)) => builder.add_leg(
                    target,
                    Flow::Sink,
                    scalar(),
                    leg("in", 0, ExternalState::Incoming, 0),
                ),
                (End::Vertex(source), End::Outgoing) => builder.add_leg(
                    source,
                    Flow::Source,
                    scalar(),
                    leg("out", 1, ExternalState::Outgoing, 0),
                ),
                _ => unreachable!("legs end at a vertex"),
            };
            edges.push((edge, source, target));
        };
        let [left, right, marker, source_side, target_side] =
            [left, right, marker, source_side, target_side].map(End::Vertex);
        if permuted_ids {
            add_edge(left, marker);
            add_edge(target_side, End::Outgoing);
            add_edge(left, right);
            add_edge(End::Incoming, source_side);
            add_edge(left, right);
        } else {
            add_edge(End::Incoming, source_side);
            add_edge(left, right);
            add_edge(left, right);
            add_edge(left, marker);
            add_edge(target_side, End::Outgoing);
        }
        let builder = builder.builder;

        let left_vertices = if global_left_right_swap {
            BTreeSet::from([End::Incoming, right])
        } else {
            BTreeSet::from([End::Incoming, left, marker])
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
        let diagram = builder.cuts(vec![cut]).build().unwrap();
        let partitions = diagram
            .cuts()
            .iter()
            .map(|cut| (cut.left.half_edges.clone(), cut.right.half_edges.clone()))
            .collect();
        let diagram = diagram
            .with_topology_threshold_partitions(partitions)
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
        let mut builder =
            FixtureBuilder::new(FeynmanDiagram::builder(Arc::clone(model), "left-right"));
        let interaction = model.vertex_rule_id("V").unwrap();
        let left = builder.add_vertex(DiagramVertex::interaction("left", interaction));
        let right = builder.add_vertex(DiagramVertex::interaction("right", interaction));
        let marker = builder.add_vertex(DiagramVertex::interaction("marker", interaction));
        let targets = match attachment {
            LeftRightAttachment::Identity => [left, left, right, right],
            LeftRightAttachment::GlobalSwap => [right, right, left, left],
            LeftRightAttachment::PartialSwap => [right, left, left, right],
        };
        let legs = [
            ("in0", 0, ExternalState::Incoming, 0),
            ("in1", 1, ExternalState::Incoming, 1),
            ("out0", 2, ExternalState::Outgoing, 0),
            ("out1", 3, ExternalState::Outgoing, 1),
        ];
        for ((name, index, state, connection), vertex) in legs.into_iter().zip(targets) {
            let flow = match state {
                ExternalState::Incoming => Flow::Sink,
                ExternalState::Outgoing => Flow::Source,
            };
            builder.add_leg(
                vertex,
                flow,
                DiagramEdge::new(particle(model, 25), false),
                leg(name, index, state, connection),
            );
        }
        builder.add_edge(left, right, DiagramEdge::new(particle(model, 25), false));
        builder.add_edge(left, marker, DiagramEdge::new(particle(model, 25), false));
        builder.builder.build().unwrap()
    }

    fn exact_options() -> GraphGroupingOptions {
        GraphGroupingOptions {
            check_canonical_numerator: true,
            ..GraphGroupingOptions::default()
        }
    }

    #[test]
    fn canonical_grouping_ignores_numeric_content_placement_in_chiral_cycles() {
        let model = Arc::new(model());
        let bridge_x = spenso::tensor_symbol!("feynkit_grouping_test::chiral_bridge_x");
        let bridge_y = spenso::tensor_symbol!("feynkit_grouping_test::chiral_bridge_y");
        let different_bridge = spenso::tensor_symbol!("feynkit_grouping_test::different_bridge");
        let spectator = test_atom("(M+N)^7");
        let closed_cycle = |indices: [usize; 7], scaled_left: bool, bridge| {
            let [a, b, c, d, x, y, mu] = indices;
            let spinor = |index: usize| Bispinor {}.to_symbolic([Atom::num(4), Atom::num(index)]);
            let lorentz = Minkowski {}.to_symbolic([Atom::num(4), Atom::num(mu)]);
            let vector = |left, right| function!(AGS.gamma, spinor(left), spinor(right), &lorentz);
            let axial = |left, right, internal| {
                vector(left, internal) * function!(AGS.gamma5, spinor(internal), spinor(right))
            };
            let scales = if scaled_left { [3, 1] } else { [1, 3] };
            // Put the number inside one chiral sum. Equality follows by moving a
            // scalar across a product, without distributing either matrix sum.
            let left = scales[0] * test_atom("M") * vector(a, b)
                + scales[0] * test_atom("N") * axial(a, b, x);
            let right = scales[1] * test_atom("M") * vector(c, d)
                + scales[1] * test_atom("N") * axial(c, d, y);
            left * right
                * function!(bridge, spinor(b), spinor(c))
                * function!(bridge_y, spinor(d), spinor(a))
                * &spectator
        };
        let grouping = Grouping::UpToScalar(GraphGroupingOptions {
            number_of_numerical_samples: 0,
            fully_numerical_substitution: true,
            ..exact_options()
        });
        for indices in [
            [701, 703, 702, 707, 704, 706, 705],
            [19, 13, 17, 11, 23, 29, 31],
        ] {
            for reversed in [false, true] {
                let mut sources: [Atom; 2] = [
                    closed_cycle([101, 102, 103, 104, 105, 106, 107], true, bridge_x),
                    closed_cycle(indices, false, bridge_x),
                ];
                if reversed {
                    sources.reverse();
                }
                let diagrams = sources
                    .iter()
                    .enumerate()
                    .map(|(index, source)| {
                        diagram_with_atoms(
                            &model,
                            &format!("cycle{index}"),
                            source.clone(),
                            Atom::one(),
                            25,
                            25,
                        )
                        .with_overall_factor(Atom::num(-1))
                    })
                    .map(DiagramPair::Shared)
                    .collect();
                let grouped = group_diagrams(
                    diagrams,
                    &model,
                    &grouping,
                    false,
                    &Run::default().threads(1),
                    false,
                    DiagramLabels::default(),
                )
                .unwrap();
                assert_eq!(
                    grouped.groups.len(),
                    1,
                    "indices={indices:?}, reversed={reversed}"
                );
                assert_eq!(grouped.groups[0].members.len(), 2);
                assert_eq!(grouped.groups[0].members[1].ratio, Atom::one());
                assert_eq!(grouped.diagrams.len(), 1);
                assert_eq!(
                    grouped.diagrams[0].numerator().unwrap().expression(),
                    &sources[0]
                );
                assert!(
                    grouped.diagrams[0]
                        .numerator()
                        .unwrap()
                        .expression()
                        .contains(&spectator)
                );
                assert_eq!(
                    evaluate_overall_factor(grouped.diagrams[0].overall_factor().as_view()),
                    Atom::num(-2),
                );
            }
        }

        // A distinct matrix factor is not an overall scalar rescaling.
        let grouped = group_diagrams(
            [bridge_x, different_bridge]
                .into_iter()
                .enumerate()
                .map(|(index, bridge)| {
                    diagram_with_atoms(
                        &model,
                        &format!("control{index}"),
                        closed_cycle([41, 42, 43, 44, 45, 46, 47], index == 0, bridge),
                        Atom::one(),
                        25,
                        25,
                    )
                })
                .map(DiagramPair::Shared)
                .collect(),
            &model,
            &grouping,
            false,
            &Run::default().threads(1),
            false,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 2);
    }

    #[test]
    fn invariant_ratios_are_scalar_without_exposing_representation_metadata() {
        let coefficient = idenso::color_cas!(2, idenso::coad!(8))
            * idenso::color_cas!(2, cof!(3)).pow(-2)
            * color_idx!(2, cof!(3));
        let scalars = BTreeSet::new();
        assert!(expression_is_scalar(&coefficient, &scalars));
        let tensor = test_atom("spenso::g(spenso::coad(8,a),spenso::coad(8,b))");
        assert!(!expression_is_scalar(&(coefficient * tensor), &scalars));
    }

    #[test]
    fn numeric_grouping_contracts_the_shared_charge_conjugation_matrix() {
        let _ = AGS.charge_conjugation;
        let model = Arc::new(model());
        let diagram = diagram(&model, "charge_conjugation", "1", 25, 25);
        let options = GraphGroupingOptions::default();
        for (expression, expected) in [
            (
                "spenso::charge_conjugation(spenso::bis(4,1),spenso::bis(4,2))*spenso::charge_conjugation(spenso::bis(4,2),spenso::bis(4,1))",
                -4,
            ),
            (
                "spenso::charge_conjugation(spenso::bis(4,1),spenso::bis(4,2))^2",
                4,
            ),
        ] {
            assert_eq!(
                evaluate_tensor_sample(&test_atom(expression), &diagram, &options, 0).unwrap(),
                Atom::num(expected),
            );
        }
    }

    #[test]
    fn sampled_lorentz_components_keep_open_color_and_zero_interfaces() {
        let model = Arc::new(model());
        let diagram = diagram(&model, "open-color-sample", "1", 25, 25);
        let options = GraphGroupingOptions::default();
        let color = feynkit_graph::NumeratorInterface::admit(
            test_atom("spenso::g(spenso::coad(8,a),spenso::coad(8,b))"),
            None,
        )
        .unwrap();
        let spin = feynkit_graph::NumeratorInterface::admit(test_atom(
            "spenso::charge_conjugation(spenso::bis(4,1),spenso::bis(4,2))*spenso::charge_conjugation(spenso::bis(4,2),spenso::bis(4,1))",
        ), None).unwrap();
        let source = color.multiply(&spin).unwrap();
        let sampled = super::evaluate_tensor_sample(&source, &diagram, &options, 0).unwrap();
        assert_eq!(sampled.rank(), 2);
        assert_eq!(sampled.structure(), color.structure());
        assert_eq!(sampled.expression(), &(Atom::num(-4) * color.expression()));
        let zero = source.with_checked_expression(Atom::Zero).unwrap();
        let sampled_zero = super::evaluate_tensor_sample(&zero, &diagram, &options, 0).unwrap();
        assert!(sampled_zero.is_zero());
        assert_eq!(sampled_zero.structure(), color.structure());
    }

    #[test]
    fn tensor_samples_preserve_symbolic_factors_and_exact_component_values() {
        let _ = AGS.charge_conjugation;
        let model = Arc::new(model());
        let diagram = diagram(&model, "factored-sample", "1", 25, 25);
        let spectator = test_atom("(UFO::M+UFO::N)^30");
        let contracted = test_atom(
            "spenso::charge_conjugation(spenso::bis(4,1),spenso::bis(4,2))*spenso::charge_conjugation(spenso::bis(4,2),spenso::bis(4,1))",
        );
        let numerator = &spectator * contracted;
        let options = GraphGroupingOptions::default();
        assert_eq!(
            evaluate_tensor_sample(&numerator, &diagram, &options, 0).unwrap(),
            Atom::num(-4) * &spectator,
        );
        let options = GraphGroupingOptions {
            fully_numerical_substitution: true,
            ..options
        };
        let m = complex_sample_value("model:UFO::M", options.numerical_sample_seed, 0);
        let n = complex_sample_value("model:UFO::N", options.numerical_sample_seed, 0);
        assert_eq!(
            evaluate_tensor_sample(&numerator, &diagram, &options, 0).unwrap(),
            Atom::num(-4) * (m + n).pow(30),
        );
    }

    #[test]
    fn finite_sample_certificate_covers_the_degree_and_rejects_unsupported_domains() {
        let left = test_atom("(x+1)*(x+2)*(x+3)");
        let right = test_atom("6+x*(11+x*(6+x))");
        assert!(PreparedNumerator::finite_samples_equal(&left, &right));
        // Six coincident probes are insufficient for this nonzero degree-six
        // polynomial. A complete certificate must also test its seventh point.
        let vanishing_probes = test_atom("x*(x-1)*(x-2)*(x-3)*(x-4)*(x-5)");
        assert!(!PreparedNumerator::finite_samples_equal(
            &vanishing_probes,
            &Atom::Zero
        ));
        assert!(!PreparedNumerator::finite_samples_equal(
            &test_atom("(x+y)^1000"),
            &Atom::Zero
        ));
        for unsupported in ["sin(x)", "(x+1)^(-1)", "(x+1)^(1/2)", "1.25`20*x"] {
            assert!(PreparedNumerator::sample_degrees(test_atom(unsupported).as_view()).is_none());
        }
        assert_eq!(left, test_atom("(x+1)*(x+2)*(x+3)"));
    }

    #[test]
    fn massive_closed_loop_flow_reversal_preserves_the_complete_nlo_inventory() {
        // The physical e+e- -> ttH NLO command's native process/options, without
        // application state. Main combines the 30 nonzero flow assignments into
        // 16 masters, including both closed-top-loop orientations in one master.
        let input: serde_json::Value =
            serde_json::from_str(include_str!("../tests/fixtures/tt_higgs_nlo_grouping.json"))
                .unwrap();
        let process: crate::Process = serde_json::from_value(input["process"].clone()).unwrap();
        let model =
            Model::from_json(include_str!("../tests/fixtures/tt_higgs_nlo_model.json")).unwrap();
        let run = Run::default().threads(1);
        let result = process
            .cross_section(crate::Loops::forward(3..=3))
            .with_orders(
                crate::Orders::new()
                    .exact("QCD", 2)
                    .exact("QED", 6)
                    .perturbative("QCD", 1),
            )
            .with_restrict(crate::Restrict::new().exclude_particles([
                1_i64, 2, 3, 4, 5, 12, 13, 14, 15, 16, 23, 24, 250, 251, 9000001, 9000002, 9000003,
                9000004,
            ]))
            .with_topology(
                crate::TopologyRules::new()
                    .with_self_energies(crate::Veto::all())
                    .with_snails(crate::Veto::massless())
                    .with_tadpoles(crate::Veto::all()),
            )
            .with_cuts(crate::Cuts::new().with_blobs(1..=1).with_spectators(0..=0))
            .generate(Arc::new(model), &run)
            .unwrap()
            .build(
                &crate::NumeratorOptions::default().with_naming(crate::Naming::prefix("GL")),
                &run,
            )
            .unwrap()
            .group(&Grouping::UpToScalar(GraphGroupingOptions::default()), &run)
            .unwrap();
        assert!(result.report().completed);
        assert_eq!(
            result.report().zero_numerator_count,
            16_u32.into(),
            "{:?}",
            result.report()
        );
        assert_eq!(result.diagrams().len(), 16, "{:?}", result.report());
        assert_eq!(
            result
                .groups()
                .iter()
                .map(|group| group.members.len())
                .sum::<usize>(),
            30
        );
    }

    #[test]
    fn finite_closed_top_loop_samples_match_under_charge_conjugation() {
        // Captured finite component coefficients for the two six-vertex closed
        // top-loop flows in e+e- -> ttH at NLO (main groups them with weight -2).
        // Their raw tensor numerators stay factorized; this certificate only
        // collects the completed exact scalar sample bodies.
        let samples: [Vec<String>; 2] = serde_json::from_str(include_str!(
            "../tests/fixtures/closed_top_loop_cc_samples.json"
        ))
        .unwrap();
        let model = Arc::new(model());
        let diagram = diagram(&model, "finite-closed-loop", "1", 25, 25);
        let options = GraphGroupingOptions::default();
        assert_eq!(samples[0].len(), 5);
        assert_eq!(samples[1].len(), 5);
        for (left, right) in samples[0].iter().zip(&samples[1]) {
            let left = evaluate_tensor_sample(&test_atom(left), &diagram, &options, 0).unwrap();
            let right = evaluate_tensor_sample(&test_atom(right), &diagram, &options, 0).unwrap();
            assert_eq!(left, right);
        }
    }

    #[test]
    fn tensor_sample_library_covers_every_finalized_external_endpoint() {
        let model = Arc::new(
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap(),
        );
        let heads = [
            feynkit_graph::symbols::epsilon(),
            feynkit_graph::symbols::epsilonbar(),
            feynkit_graph::symbols::u(),
            feynkit_graph::symbols::ubar(),
            feynkit_graph::symbols::v(),
            feynkit_graph::symbols::vbar(),
        ];
        let spectator = test_atom("(UFO::MT+UFO::MH)^11");
        for pdg in [22, 6, -6] {
            for paired in [false, true] {
                for state in [ExternalState::Incoming, ExternalState::Outgoing] {
                    for reverse_storage in [false, true] {
                        for connection in [0, 3] {
                            let mut builder =
                                FeynmanDiagram::builder(model.clone(), "sample-endpoints");
                            let vertices = ["left", "right"].map(|name| {
                                builder.add_vertex(DiagramVertex {
                                    name: name.into(),
                                    interaction: None,
                                    numerator: SymbolicTensor::infer(Atom::one()).unwrap(),
                                })
                            });
                            let mut edge = DiagramEdge::new(particle(&model, pdg), pdg != 22);
                            if paired {
                                edge.external = Some(ExternalLeg {
                                    name: "paired".into(),
                                    index: 0,
                                    state,
                                    connection,
                                });
                                let (source, target) = if reverse_storage {
                                    (vertices[1], vertices[0])
                                } else {
                                    (vertices[0], vertices[1])
                                };
                                builder.add_edge(source, target, edge).unwrap();
                            } else {
                                // A generated leg stored against its process flow is
                                // conjugated when finalization orients it.
                                let flow = match state {
                                    ExternalState::Incoming => Flow::Sink,
                                    ExternalState::Outgoing => Flow::Source,
                                };
                                let flow = if reverse_storage { -flow } else { flow };
                                builder
                                    .add_generation_external(
                                        vertices[0],
                                        flow,
                                        edge,
                                        VertexSlot(0),
                                        ExternalLeg {
                                            name: "external".into(),
                                            index: 0,
                                            state,
                                            connection,
                                        },
                                    )
                                    .unwrap();
                                builder
                                    .add_generation_external(
                                        vertices[1],
                                        Flow::Source,
                                        DiagramEdge::new(particle(&model, 25), false),
                                        VertexSlot(0),
                                        ExternalLeg {
                                            name: "spectator".into(),
                                            index: 1,
                                            state: ExternalState::Outgoing,
                                            connection: connection + 17,
                                        },
                                    )
                                    .unwrap();
                            }
                            builder
                                .add_edge(
                                    vertices[0],
                                    vertices[1],
                                    DiagramEdge::new(particle(&model, 25), false),
                                )
                                .unwrap();
                            let diagram = builder.build().unwrap();
                            if paired && state == ExternalState::Outgoing {
                                assert!(matches!(
                                    diagram.validate(),
                                    Err(DiagramError::Invariant {
                                        operation: "validating sewn external flow",
                                        ..
                                    })
                                ));
                                continue;
                            }
                            let diagram = if paired {
                                let (left, right) = diagram.half_edges().partition(|half_edge| {
                                    let node = diagram
                                        .underlying()
                                        .node_id(diagram.half_edge_id(*half_edge).unwrap());
                                    VertexId(node.0) == vertices[0]
                                });
                                diagram.with_cut_partitions(vec![(left, right)]).unwrap()
                            } else {
                                diagram
                            };
                            diagram.validate().unwrap();
                            // This independent producer owns endpoint semantics.
                            let projector = diagram.external_state_projector().unwrap();
                            let projector = feynkit_graph::NumeratorInterface::export(&projector).0;
                            let expected_heads = heads
                                .into_iter()
                                .filter(|head| projector.contains_symbol(*head))
                                .collect::<Vec<_>>();
                            assert_eq!(expected_heads.len(), if paired { 2 } else { 1 });
                            for symmetric_polarizations in [false, true] {
                                let options = GraphGroupingOptions {
                                    symmetric_polarizations,
                                    ..GraphGroupingOptions::default()
                                };
                                let library = tensor_sample_library(&diagram, &options, 0).unwrap();
                                let mut values = Vec::new();
                                for head in heads {
                                    let representation = if heads[..2].contains(&head) {
                                        Minkowski {}.new_rep(4).cast::<LibraryRep>()
                                    } else {
                                        Bispinor {}.new_rep(4).cast::<LibraryRep>()
                                    };
                                    let key = ExplicitKey::from_iter(
                                        [representation],
                                        head,
                                        Some(vec![
                                            symbol!("FeynKit::ExternalMomentum")
                                                .call(connection as i64),
                                        ]),
                                    );
                                    let tensor = library.get_storage(key.canonical());
                                    assert_eq!(
                                        tensor.is_ok(),
                                        expected_heads.contains(&head),
                                        "pdg={pdg}, paired={paired}, state={state:?}, reverse={reverse_storage}, head={head}"
                                    );
                                    if let Ok(tensor) = tensor {
                                        let actual = tensor
                                            .iter_flat()
                                            .map(|(_, value)| value.to_owned())
                                            .collect::<Vec<Atom>>();
                                        let label = if symmetric_polarizations {
                                            format!("polarization:{connection}")
                                        } else {
                                            format!("polarization:{head}:{connection}")
                                        };
                                        let expected = (0..4)
                                            .map(|component| {
                                                complex_sample_value(
                                                    &format!("{label}:{component}"),
                                                    options.numerical_sample_seed,
                                                    0,
                                                )
                                            })
                                            .collect::<Vec<_>>();
                                        assert_eq!(actual, expected);
                                        values.push(actual);
                                    }
                                    let foreign = ExplicitKey::from_iter(
                                        [representation],
                                        head,
                                        Some(vec![
                                            symbol!("FeynKit::ExternalMomentum")
                                                .call((connection + 1) as i64),
                                        ]),
                                    );
                                    assert!(library.get(foreign.canonical()).is_err());
                                }
                                if paired {
                                    if symmetric_polarizations {
                                        assert_eq!(values[0], values[1]);
                                    } else {
                                        assert_ne!(values[0], values[1]);
                                        assert_ne!(
                                            values[0],
                                            values[1]
                                                .iter()
                                                .map(AtomCore::conj)
                                                .collect::<Vec<_>>()
                                        );
                                    }
                                    let mut ports = Vec::new();
                                    projector.visitor(&mut |atom| {
                                        if let AtomView::Fun(function) = atom
                                            && expected_heads.contains(&function.get_symbol())
                                        {
                                            ports.push(function.get(1).to_owned());
                                            false
                                        } else {
                                            true
                                        }
                                    });
                                    assert_eq!(ports.len(), 2);
                                    let numerator =
                                        &spectator * &projector * spenso::g!(&ports[0], &ports[1]);
                                    let routed = normalize_momentum_routing(&diagram, numerator);
                                    let scalar = Atom::add_many((0..4).map(|i| {
                                        let product = &values[0][i] * &values[1][i];
                                        if pdg == 22 && i != 0 {
                                            -product
                                        } else {
                                            product
                                        }
                                    }));
                                    assert_eq!(
                                        evaluate_tensor_sample(&routed, &diagram, &options, 0)
                                            .unwrap(),
                                        &scalar * &spectator,
                                    );
                                    let numeric_options = GraphGroupingOptions {
                                        fully_numerical_substitution: true,
                                        ..options
                                    };
                                    let mt = complex_sample_value(
                                        "model:UFO::MT",
                                        options.numerical_sample_seed,
                                        0,
                                    );
                                    let mh = complex_sample_value(
                                        "model:UFO::MH",
                                        options.numerical_sample_seed,
                                        0,
                                    );
                                    assert_eq!(
                                        evaluate_tensor_sample(
                                            &routed,
                                            &diagram,
                                            &numeric_options,
                                            0
                                        )
                                        .unwrap(),
                                        scalar * (mt + mh).pow(11),
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn charge_conjugate_sewn_self_energy_samples_bind_both_endpoints() {
        let _ = AGS.gamma;
        let model = Arc::new(
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap(),
        );
        // Actual pre-UV native inputs for the two massive-quark self-energy
        // supergraphs in gamma* -> t tbar. Preserve their matrix products.
        let sources = [
            test_atom(include_str!(
                "../tests/fixtures/sewn_top_self_energy_cc_left.atom"
            )),
            test_atom(include_str!(
                "../tests/fixtures/sewn_top_self_energy_cc_right.atom"
            )),
        ];
        // C gamma^T C^-1=-gamma, so a closed spinor circuit is unchanged by
        // transposing every matrix and inserting these signs. This exact
        // certificate is independent of the numerical grouping decision.
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        for source in &sources {
            let tensor =
                SymbolicTensor::infer(cooking.try_cook(source.as_view()).unwrap()).unwrap();
            assert_eq!(
                tensor
                    .structure()
                    .canonical()
                    .external_structure_iter()
                    .count(),
                0
            );
        }
        assert_eq!(conjugate_spinor_circuits(&sources[0]), sources[1]);
        assert_eq!(
            conjugate_spinor_circuits(&conjugate_spinor_circuits(&sources[0])),
            sources[0]
        );

        let diagrams = sources
            .into_iter()
            .enumerate()
            .map(|(index, numerator)| {
                let mut builder =
                    FeynmanDiagram::builder(model.clone(), format!("self_energy_{index}"));
                let vertices = ["qcd0", "qcd1", "photon0", "photon1"].map(|name| {
                    builder.add_vertex(DiagramVertex {
                        name: name.into(),
                        interaction: None,
                        numerator: SymbolicTensor::infer(Atom::one()).unwrap(),
                    })
                });
                let mut external = DiagramEdge::new(particle(&model, 22), false);
                external.external = Some(ExternalLeg {
                    name: "ext0".into(),
                    index: 0,
                    state: ExternalState::Incoming,
                    connection: 0,
                });
                builder
                    .add_edge(vertices[3], vertices[2], external)
                    .unwrap();
                let sign = if index == 0 { -1 } else { 1 };
                for (source, target, pdg) in [
                    (0, 1, sign * 6),
                    (0, 1, 21),
                    (0, 3, -sign * 6),
                    (1, 2, sign * 6),
                    (2, 3, sign * 6),
                ] {
                    builder
                        .add_edge(
                            vertices[source],
                            vertices[target],
                            DiagramEdge::new(particle(&model, pdg), pdg != 21),
                        )
                        .unwrap();
                }
                let diagram = builder.build().unwrap();
                let (left, right) = diagram.half_edges().partition(|half_edge| {
                    let node = diagram
                        .underlying()
                        .node_id(diagram.half_edge_id(*half_edge).unwrap());
                    VertexId(node.0) == vertices[2]
                });
                let diagram = diagram.with_cut_partitions(vec![(left, right)]).unwrap();
                // Replay the recorded grouping input without claiming its
                // replacement fragment is a local Feynman-rule numerator.
                let diagram = diagram
                    .with_numerator(
                        feynkit_graph::NumeratorInterface::admit(numerator, None).unwrap(),
                    )
                    .unwrap();
                diagram.validate().unwrap();
                diagram
            })
            .collect::<Vec<_>>();
        let options = GraphGroupingOptions::default();
        assert!(scalar_names(&model).contains("MT"));
        assert_eq!(
            diagrams[0].loop_momentum_basis(),
            diagrams[1].loop_momentum_basis()
        );
        assert_eq!(
            canonical_topology(&diagrams[0], &model, &options, true)
                .unwrap()
                .key,
            canonical_topology(&diagrams[1], &model, &options, true)
                .unwrap()
                .key,
        );
        let mut sampled_words = Vec::new();
        for diagram in &diagrams {
            let samples =
                numerical_tensor_samples(diagram.numerator().unwrap(), diagram, &options).unwrap();
            for value in &samples {
                assert!(
                    !value
                        .expression()
                        .contains_symbol(feynkit_graph::symbols::epsilon())
                );
                assert!(
                    !value
                        .expression()
                        .contains_symbol(feynkit_graph::symbols::epsilonbar())
                );
                assert!(value.expression().contains_symbol(symbol!("UFO::MT")));
            }
            sampled_words.push(samples);
        }
        for (sample, (left, right)) in sampled_words[0].iter().zip(&sampled_words[1]).enumerate() {
            // The source words obey the exact charge-conjugation identity
            // above. These finite component probes also check its numerical
            // execution without distributing the recorded numerator.
            for mass in 0..=4 {
                let left = left
                    .with_checked_expression(
                        left.expression()
                            .replace(symbol!("UFO::MT"))
                            .with(Atom::num(mass)),
                    )
                    .unwrap();
                let right = right
                    .with_checked_expression(
                        right
                            .expression()
                            .replace(symbol!("UFO::MT"))
                            .with(Atom::num(mass)),
                    )
                    .unwrap();
                assert!(
                    super::expressions_equal(&left, &right).unwrap(),
                    "sample={sample}, mass={mass}: left={left}; right={right}"
                );
            }
            assert!(
                super::expressions_equal(left, right).unwrap(),
                "symbolic sample={sample}: left={left}; right={right}"
            );
        }
    }

    #[test]
    fn factor_preserving_comparison_retains_unproved_identities_and_checks_tensor_scopes() {
        let left = test_atom("x*(y+z)");
        let right = test_atom("x*y+x*z");
        assert!(expressions_equal(&left, &right).unwrap());
        assert!(expressions_equal(&test_atom("x+y"), &test_atom("y+x")).unwrap());
        assert!(!expressions_equal(&test_atom("(x+y)^2"), &test_atom("x^2+2*x*y+y^2"),).unwrap());
        assert!(
            expressions_equal(
                &test_atom("v(spenso::mink(4,1))^3"),
                &test_atom("v(spenso::mink(4,1))^3"),
            )
            .is_err()
        );
    }

    #[test]
    fn grouping_progress_counts_zero_and_retained_inputs() {
        let model = Arc::new(model());
        let updates = Arc::new(std::sync::Mutex::new(Vec::new()));
        let captured = updates.clone();
        let options = Run::default().progress(move |snapshot| {
            captured.lock().unwrap().push(snapshot);
            crate::GenerationControl::Continue
        });
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x+y", 25, 25),
                diagram(&model, "g1", "y+x", 25, 25),
                diagram(&model, "g2", "0", 25, 25),
            ]
            .into_iter()
            .map(DiagramPair::Shared)
            .collect(),
            &model,
            &Grouping::Identical(exact_options()),
            false,
            &options,
            false,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.zero_numerator_count, 1);
        let updates = updates.lock().unwrap();
        for (stage, total) in [
            (crate::Stage::GroupingPreparation, 3),
            (crate::Stage::GroupingSamples, 2),
            (crate::Stage::GroupingComparison, 2),
        ] {
            let counts = updates
                .iter()
                .filter(|snapshot| snapshot.stage == stage)
                .map(|snapshot| {
                    (
                        usize::try_from(&snapshot.completed).unwrap(),
                        snapshot.total.as_ref().map(|n| usize::try_from(n).unwrap()),
                    )
                })
                .collect::<Vec<_>>();
            assert_eq!(
                counts,
                (0..=total).map(|n| (n, Some(total))).collect::<Vec<_>>()
            );
        }
        assert_eq!(updates.last().unwrap().stage, crate::Stage::Grouping);
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn transient_cancellation_never_compares_unprepared_numerators() {
        let model = Arc::new(model());
        let inputs = vec![
            diagram(&model, "first", "2", 25, 25),
            diagram(&model, "second", "1", 25, 25),
        ];
        let expected = inputs
            .iter()
            .map(|diagram| diagram.to_json().unwrap())
            .collect::<Vec<_>>();
        let signal = Arc::new(AtomicBool::new(false));
        let arm = signal.clone();
        let poll = signal.clone();
        let stages = Arc::new(Mutex::new(Vec::new()));
        let observed = stages.clone();
        let options = Run::default()
            .progress(move |update| {
                observed.lock().unwrap().push(update.stage);
                if update.stage == crate::Stage::ZeroColor && update.completed == 0usize.into() {
                    arm.store(true, Ordering::Release);
                }
                crate::GenerationControl::Continue
            })
            .cancellation_check(move || poll.swap(false, Ordering::AcqRel));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let grouped = pool
            .install(|| {
                group_diagrams(
                    inputs.into_iter().map(DiagramPair::Shared).collect(),
                    &model,
                    &Grouping::Identical(exact_options()),
                    false,
                    &options,
                    true,
                    DiagramLabels::default(),
                )
            })
            .unwrap();

        // Skipping the first numerator leaves an internal placeholder of one.
        // A cleared signal must not let it compare equal to the second input.
        assert!(!signal.load(Ordering::Acquire));
        assert!(options.cancellation_requested());
        assert_eq!(grouped.zero_numerator_count, 0);
        assert_eq!(
            grouped
                .diagrams
                .iter()
                .map(|diagram| diagram.to_json().unwrap())
                .collect::<Vec<_>>(),
            expected
        );
        assert_eq!(grouped.groups.len(), 2);
        for (source, group) in grouped.groups.iter().enumerate() {
            assert_eq!(group.master, source);
            assert_eq!(group.members.len(), 1);
            assert_eq!(group.members[0].source_diagram, source);
            assert_eq!(group.members[0].diagram, source);
        }
        assert!(
            stages
                .lock()
                .unwrap()
                .iter()
                .all(|&stage| stage == crate::Stage::ZeroColor)
        );
    }

    #[test]
    fn shared_finalized_masters_move_through_grouping() {
        let model = Arc::new(model());
        for grouping in [Grouping::None, Grouping::Identical(exact_options())] {
            let first = diagram(&model, "first", "1", 25, 25);
            let allocation = first.vertices().next().unwrap().1 as *const DiagramVertex;
            let second = diagram(&model, "second", "1", 25, 25);
            let grouped = group_diagrams(
                vec![DiagramPair::Shared(first), DiagramPair::Shared(second)],
                &model,
                &grouping,
                false,
                &Run::default(),
                false,
                DiagramLabels::default(),
            )
            .unwrap();
            assert_eq!(
                grouped.diagrams[0].vertices().next().unwrap().1 as *const DiagramVertex,
                allocation,
                "the retained master must keep its owned graph allocation",
            );
            assert!(
                grouped
                    .groups
                    .iter()
                    .flat_map(|group| &group.members)
                    .all(|member| member.source_identity.is_none())
            );
            assert_eq!(
                grouped.diagrams[0].numerator().unwrap().expression(),
                &Atom::one()
            );
        }
    }

    #[test]
    fn cancelled_grouping_moves_the_distinct_physical_owner() {
        let model = Arc::new(model());
        let comparison = diagram(&model, "source", "0", 25, 25);
        let representative = diagram(&model, "source", "x+y", 25, 25);
        let expected = representative.to_json().unwrap();
        let allocation = representative.vertices().next().unwrap().1 as *const DiagramVertex;
        let grouped = group_diagrams(
            vec![DiagramPair::Distinct {
                comparison,
                representative: Box::new(representative),
                reversed_edges: BTreeSet::new(),
            }],
            &model,
            &Grouping::Identical(exact_options()),
            false,
            &Run::default().cancellation_check(|| true),
            false,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.zero_numerator_count, 0);
        assert_eq!(grouped.diagrams.len(), 1);
        assert_eq!(grouped.diagrams[0].to_json().unwrap(), expected);
        assert_eq!(
            grouped.diagrams[0].vertices().next().unwrap().1 as *const DiagramVertex,
            allocation
        );
    }

    #[test]
    fn groups_identical_numerators() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x+y", 25, 25),
                diagram(&model, "g1", "y+x", 25, 25),
            ]
            .into_iter()
            .map(DiagramPair::Shared)
            .collect(),
            &model,
            &Grouping::Identical(exact_options()),
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
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
            vec![external_first, internal_first]
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
            &model,
            &Grouping::Identical(GraphGroupingOptions {
                check_canonical_numerator: false,
                number_of_numerical_samples: 2,
                ..GraphGroupingOptions::default()
            }),
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
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
                diagram(&model, "g0", "x+y", 25, 25).with_overall_factor(Atom::num(2)),
                diagram(&model, "g1", "-x-y", 25, 25),
            ]
            .into_iter()
            .map(DiagramPair::Shared)
            .collect(),
            &model,
            &Grouping::UpToSign(exact_options()),
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
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
            ]
            .into_iter()
            .map(DiagramPair::Shared)
            .collect(),
            &model,
            &Grouping::UpToScalar(exact_options()),
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members[1].ratio, test_atom("2"));
    }

    // Frozen DIS19 local factors and signed routes; only legacy particle IDs
    // were translated to names. The restricted SM keeps the source mass classes.
    fn sampled_quotient_physical_pair() -> (Arc<Model>, Vec<FeynmanDiagram>) {
        let mut model = Model::standard_model();
        let card = ParameterCard::from_path(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../assets/models/json/sm/restrict_default.json"
        ))
        .unwrap();
        model.apply_parameter_card(&card).unwrap();
        model.simplify_zero_values().unwrap();
        let model = Arc::new(model);
        let diagrams = [
            include_str!("../tests/fixtures/sampled_quotient/GL29.dot"),
            include_str!("../tests/fixtures/sampled_quotient/GL68.dot"),
        ]
        .map(|dot| FeynmanDiagram::from_dot(Arc::clone(&model), dot).unwrap())
        .to_vec();
        for diagram in &diagrams {
            for (_, _, edge) in diagram.edges() {
                let mass = model.particle_by_id(edge.particle).unwrap().mass;
                assert_eq!(model.parameter_by_id(mass).unwrap().name, "ZERO");
            }
            assert_eq!(
                evaluate_overall_factor(diagram.overall_factor().as_view()),
                -Atom::one()
            );
        }
        let [master, candidate] = diagrams.as_slice() else {
            unreachable!()
        };
        assert_eq!(
            master.loop_momentum_basis(),
            candidate.loop_momentum_basis()
        );
        assert_eq!(master.projector(), candidate.projector());
        assert_eq!(
            master.numerator_prefactor(),
            candidate.numerator_prefactor()
        );
        for ((id, vertex), (other_id, other)) in master.vertices().zip(candidate.vertices()) {
            assert_eq!(id, other_id);
            let ratio = Atom::num(if id.0 == 0 { -2_i64 } else { 1 });
            assert_eq!(
                model
                    .expand_couplings(other.numerator.expression())
                    .collect_factors(),
                (ratio * model.expand_couplings(vertex.numerator.expression())).collect_factors(),
                "exact local coupling/factor witness at vertex {id:?}",
            );
        }
        for ((id, endpoints, edge), (other_id, other_endpoints, other)) in
            master.edges().zip(candidate.edges())
        {
            assert_eq!((id, endpoints), (other_id, other_endpoints));
            assert_eq!(edge.numerator, other.numerator);
        }
        (model, diagrams)
    }

    #[test]
    fn sampled_quotient_groups_physical_local_factors() {
        const REGISTRATION: &str = "FEYNKIT_QUOTIENT_REGISTRATION";
        if let Ok(case) = env::var(REGISTRATION) {
            let (order, threads) = case.split_once(':').unwrap();
            let threads = threads.parse::<usize>().unwrap();
            assert!(matches!(order, "forward" | "reverse"));
            assert!([1, 8, 32].contains(&threads));
            // #209 preserves valid native Named indices. Deliberately construct
            // public reversible aliases, and use them in the actual physical inputs.
            let cooking = CookSettings::indices()
                .with_source_filter(CookSourceFilter::AnyFunction)
                .with_mode(CookMode::ReversibleEncoding);
            let mut labels = (0..18).collect::<Vec<i64>>();
            if order == "reverse" {
                labels.reverse();
            }
            let mut registered = BTreeMap::new();
            let mut replacements = Vec::new();
            for &label in &labels {
                let payload = function!(symbols::hedge_index(), label, 1);
                let AtomView::Fun(function) = payload.as_view() else {
                    unreachable!()
                };
                let cooked = cooking.cook_function(function).unwrap();
                let AtomView::Var(variable) = cooked.as_view() else {
                    panic!("explicit codec did not create an index symbol")
                };
                registered.insert(variable.get_symbol(), label);
                replacements.push(Replacement::new(payload, cooked));
            }
            assert_eq!(registered.values().copied().collect::<Vec<_>>(), labels);
            let (model, diagrams) = sampled_quotient_physical_pair();
            let encode = |value: &SymbolicTensor<PartialStructure>| {
                NumeratorInterface::translate::<DiagramError>(value, |atom| {
                    Ok(atom.replace_multiple(&replacements))
                })
                .unwrap()
            };
            let diagrams = diagrams
                .into_iter()
                .map(|diagram| {
                    let vertices = diagram
                        .vertices()
                        .map(|(id, vertex)| (id, encode(&vertex.numerator)))
                        .collect::<Vec<_>>();
                    let edges = diagram
                        .edges()
                        .map(|(id, _, edge)| (id, encode(&edge.numerator)))
                        .collect::<Vec<_>>();
                    let prefactor = encode(diagram.numerator_prefactor());
                    let projector = encode(diagram.projector());
                    let encoded = diagram
                        .with_numerator_prefactor(prefactor)
                        .with_projector(projector)
                        .with_numerator_fragments(vertices, edges)
                        .unwrap();
                    let mut used = BTreeSet::new();
                    for value in [encoded.numerator().unwrap(), encoded.projector()] {
                        value.expression().visitor(&mut |part| {
                            if let AtomView::Var(variable) = part {
                                let symbol = variable.get_symbol();
                                if registered.contains_key(&symbol) {
                                    used.insert(symbol);
                                }
                            }
                            true
                        });
                    }
                    assert_eq!(
                        used.len(),
                        labels.len(),
                        "all18 actual hedge aliases must occur in the physical input"
                    );
                    encoded
                })
                .collect::<Vec<_>>();
            let master_numerator = diagrams[0].numerator().unwrap().clone();
            let options = GraphGroupingOptions {
                check_canonical_numerator: false,
                ..Default::default()
            };
            let observed = Arc::new(AtomicUsize::new(0));
            let sink = Arc::clone(&observed);
            let caller = thread::current().id();
            let run = Run::default().threads(threads).progress(move |progress| {
                if progress.stage == Stage::GroupingPreparation
                    && progress.completed > 0usize.into()
                {
                    if threads == 1 {
                        assert_eq!(thread::current().id(), caller);
                        assert!(rayon::current_thread_index().is_none());
                    } else {
                        assert_eq!(rayon::current_num_threads(), threads);
                    }
                    sink.fetch_add(1, Ordering::Relaxed);
                }
                GenerationControl::Continue
            });
            let grouped = run
                .in_pool(|| {
                    group_diagrams(
                        diagrams.into_iter().map(DiagramPair::Shared).collect(),
                        &model,
                        &Grouping::UpToScalar(options),
                        true,
                        &run,
                        false,
                        DiagramLabels::default(),
                    )
                })
                .unwrap()
                .unwrap();
            assert!(
                observed.load(Ordering::Relaxed) > 0,
                "physical witness did not observe its actual preparation pool"
            );
            assert_eq!(grouped.groups.len(), 1, "threads={threads}");
            assert_eq!(grouped.groups[0].members.len(), 2);
            assert_eq!(grouped.groups[0].members[0].ratio, Atom::one());
            assert_eq!(grouped.groups[0].members[1].ratio, Atom::num(-2));
            assert_eq!(
                evaluate_overall_factor(grouped.diagrams[0].overall_factor().as_view()),
                Atom::one()
            );
            assert_eq!(grouped.diagrams[0].numerator().unwrap(), &master_numerator);
            let members = grouped.groups[0]
                .members
                .iter()
                .map(|member| {
                    (
                        member.source_diagram,
                        member.ratio.to_canonical_string(),
                        member.overall_factor.to_canonical_string(),
                    )
                })
                .collect::<Vec<_>>();
            let graph: Value =
                serde_json::from_str(&grouped.diagrams[0].to_json().unwrap()).unwrap();
            println!(
                "FEYNKIT_QUOTIENT_RESULT {}",
                serde_json::to_string(&(graph, members)).unwrap()
            );
            return;
        }
        let mut reference = None;
        for order in ["forward", "reverse"] {
            for threads in [1, 8, 32] {
                let case = format!("{order}:{threads}");
                let output = Command::new(env::current_exe().unwrap())
                    .args([
                        "--exact",
                        "grouping::tests::sampled_quotient_groups_physical_local_factors",
                        "--nocapture",
                    ])
                    .env(REGISTRATION, &case)
                    .output()
                    .unwrap();
                assert!(
                    output.status.success()
                        && String::from_utf8_lossy(&output.stdout)
                            .contains("1 passed; 0 failed; 0 ignored"),
                    "{case}: {}\n{}\n{}",
                    output.status,
                    String::from_utf8_lossy(&output.stdout),
                    String::from_utf8_lossy(&output.stderr),
                );
                let text = String::from_utf8(output.stdout).unwrap();
                let snapshot = text
                    .lines()
                    .find(|line| line.starts_with("FEYNKIT_QUOTIENT_RESULT "))
                    .expect(
                        "fresh child did not export its physical graph and exact group members",
                    );
                if let Some(reference) = &reference {
                    assert_eq!(snapshot, reference, "{case}: complete graph/group payload");
                } else {
                    reference = Some(snapshot.to_owned());
                }
            }
        }
    }

    #[test]
    fn sampled_quotient_requires_scalar_interfaces_and_all_samples() {
        #[derive(Deserialize)]
        struct Sample {
            expression: String,
            interface: Option<String>,
        }
        let frozen: [Vec<Sample>; 2] = serde_json::from_str(include_str!(
            "../tests/fixtures/sampled_quotient/prepared_samples.json"
        ))
        .unwrap();
        let [mut master, mut candidate] = frozen.map(|samples| PreparedNumerator {
            exact: SymbolicTensor::infer(Atom::one()).unwrap(),
            sample_source: SymbolicTensor::infer(Atom::one()).unwrap(),
            samples: samples
                .into_iter()
                .map(|sample| {
                    let interface = sample
                        .interface
                        .map(|ports| NumeratorInterface::parse(&ports).unwrap());
                    NumeratorInterface::admit(test_atom(&sample.expression), interface).unwrap()
                })
                .collect(),
        });
        let scalars = scalar_names(&Model::standard_model());
        let options = GraphGroupingOptions {
            check_canonical_numerator: false,
            ..Default::default()
        };
        let scalar_ratio = |candidate: &PreparedNumerator, master: &PreparedNumerator| {
            compare(
                candidate,
                master,
                ComparisonMode::UpToScalar,
                &options,
                &scalars,
            )
            .unwrap()
        };
        assert_eq!(master.samples.len(), options.number_of_numerical_samples);
        assert_eq!(scalar_ratio(&candidate, &master), Some(Atom::num(-2)));
        assert!(
            compare(
                &candidate,
                &master,
                ComparisonMode::UpToSign,
                &options,
                &scalars
            )
            .unwrap()
            .is_none()
        );
        let original = candidate.samples.clone();
        candidate.samples = master
            .samples
            .iter()
            .map(|value| value.with_checked_expression(-value.expression()).unwrap())
            .collect();
        assert_eq!(
            compare(
                &candidate,
                &master,
                ComparisonMode::UpToSign,
                &options,
                &scalars
            )
            .unwrap(),
            Some(-Atom::one())
        );

        // All earlier samples agree: a disagreement only in the final sample
        // must still reject the proposed common ratio.
        candidate.samples = original.clone();
        candidate.samples[4] = candidate.samples[4]
            .with_checked_expression(3 * candidate.samples[4].expression())
            .unwrap();
        assert!(scalar_ratio(&candidate, &master).is_none());

        // Perturb one of the two actual color-basis coefficients independently.
        candidate.samples = original.clone();
        let sectors = candidate.samples[0]
            .coefficient_list(TensorCollectFilter::Reps([
                ColorAdjoint {}.into(),
                ColorFundamental {}.into(),
                ColorSextet {}.into(),
            ]))
            .unwrap();
        assert_eq!(sectors.len(), 2);
        candidate.samples[0] = candidate.samples[0]
            .with_checked_expression(candidate.samples[0].expression() + sectors[0].0.expression())
            .unwrap();
        assert!(scalar_ratio(&candidate, &master).is_none());

        candidate.samples = original.clone();
        // Equal zero expressions must retain distinct declared interfaces;
        // otherwise the finite-sample equality shortcut would accept them.
        candidate.samples[4] = NumeratorInterface::admit(
            test_atom("spenso::g(spenso::cof(3,other_a),spenso::dind(spenso::cof(3,other_b)))"),
            None,
        )
        .unwrap()
        .with_algebra_result(Atom::Zero)
        .unwrap();
        master.samples[4] = master.samples[4].with_algebra_result(Atom::Zero).unwrap();
        assert!(candidate.samples[4].expression().is_zero());
        assert!(master.samples[4].expression().is_zero());
        assert!(!InterfaceInference::additive_interfaces_match(
            candidate.samples[4].structure(),
            master.samples[4].structure(),
        ));
        assert!(scalar_ratio(&candidate, &master).is_none());

        // A zero master sample cannot define a nonzero candidate/master ratio.
        candidate.samples = original;
        assert!(scalar_ratio(&candidate, &master).is_none());
    }

    #[test]
    fn sampled_quotient_rejects_undefined_ratios() {
        let scalars = scalar_names(&model());
        for (candidate, master) in [
            (Atom::one(), Atom::Zero),
            (Atom::one() / Atom::Zero, Atom::one()),
            (Atom::Zero / Atom::Zero, Atom::one()),
        ] {
            assert!(
                compare_atoms(&candidate, &master, ComparisonMode::UpToScalar, &scalars)
                    .unwrap()
                    .is_none()
            );
        }
    }

    #[test]
    fn sampled_quotient_does_not_group_nonfinite_operands() {
        let model = Arc::new(model());
        for canonical in [false, true] {
            let options = GraphGroupingOptions {
                check_canonical_numerator: canonical,
                ..Default::default()
            };
            for nonfinite in [Atom::one() / Atom::Zero, Atom::Zero / Atom::Zero] {
                for (master, candidate) in [
                    (nonfinite.clone(), nonfinite.clone()),
                    (nonfinite.clone(), Atom::one()),
                    (Atom::one(), nonfinite.clone()),
                ] {
                    for grouping in [
                        Grouping::Identical(options.clone()),
                        Grouping::UpToSign(options.clone()),
                        Grouping::UpToScalar(options.clone()),
                    ] {
                        let diagrams = [master.clone(), candidate.clone()]
                            .into_iter()
                            .enumerate()
                            .map(|(index, value)| {
                                DiagramPair::Shared(diagram_with_atoms(
                                    &model,
                                    &format!("nonfinite-{index}"),
                                    value,
                                    Atom::one(),
                                    25,
                                    25,
                                ))
                            })
                            .collect();
                        let grouped = group_diagrams(
                            diagrams,
                            &model,
                            &grouping,
                            false,
                            &Run::default(),
                            false,
                            DiagramLabels::default(),
                        )
                        .unwrap();
                        assert_eq!(grouped.zero_numerator_count, 0);
                        assert_eq!(grouped.groups.len(), 2, "{candidate}/{nonfinite}");
                    }
                }
            }
        }
    }

    #[test]
    fn sampled_quotient_checks_nonfinite_operands_in_every_sample() {
        for canonical in [false, true] {
            let options = GraphGroupingOptions {
                check_canonical_numerator: canonical,
                ..Default::default()
            };
            let scalars = scalar_names(&model());
            let prepared = |expression: Atom| {
                let value = NumeratorInterface::admit(expression, None).unwrap();
                PreparedNumerator {
                    exact: value.clone(),
                    sample_source: value.clone(),
                    samples: vec![value; options.number_of_numerical_samples],
                }
            };
            for mode in [
                ComparisonMode::Identical,
                ComparisonMode::UpToSign,
                ComparisonMode::UpToScalar,
            ] {
                // Finite equality, including zero/zero, remains a valid certificate.
                for finite in [Atom::one(), Atom::Zero] {
                    assert_eq!(
                        compare(
                            &prepared(finite.clone()),
                            &prepared(finite),
                            mode,
                            &options,
                            &scalars,
                        )
                        .unwrap(),
                        Some(Atom::one())
                    );
                }
                for nonfinite in [Atom::one() / Atom::Zero, Atom::Zero / Atom::Zero] {
                    for (candidate, master) in [
                        (nonfinite.clone(), nonfinite.clone()),
                        (Atom::one(), nonfinite.clone()),
                        (nonfinite, Atom::one()),
                    ] {
                        for sample in 0..options.number_of_numerical_samples {
                            let mut candidate_samples = prepared(Atom::one());
                            let mut master_samples = prepared(Atom::one());
                            candidate_samples.samples[sample] =
                                NumeratorInterface::admit(candidate.clone(), None).unwrap();
                            master_samples.samples[sample] =
                                NumeratorInterface::admit(master.clone(), None).unwrap();
                            assert!(
                                compare(
                                    &candidate_samples,
                                    &master_samples,
                                    mode,
                                    &options,
                                    &scalars,
                                )
                                .unwrap()
                                .is_none(),
                                "sample {sample}: {candidate}/{master}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn canonical_quotient_rejects_nonfinite_operands_before_equality() {
        let scalars = scalar_names(&model());
        for mode in [
            ComparisonMode::Identical,
            ComparisonMode::UpToSign,
            ComparisonMode::UpToScalar,
        ] {
            for finite in [Atom::Zero, Atom::one()] {
                assert_eq!(
                    compare_atoms(&finite, &finite, mode, &scalars).unwrap(),
                    Some(Atom::one())
                );
            }
            for nonfinite in [Atom::one() / Atom::Zero, Atom::Zero / Atom::Zero] {
                for (candidate, master) in [
                    (nonfinite.clone(), nonfinite.clone()),
                    (Atom::one(), nonfinite.clone()),
                    (nonfinite, Atom::one()),
                ] {
                    assert!(
                        compare_atoms(&candidate, &master, mode, &scalars)
                            .unwrap()
                            .is_none()
                    );
                }
            }
        }
    }

    #[test]
    fn canonical_quotient_rejects_nonfinite_concrete_dimensions() {
        let nc = Atom::var(CS.nc);
        let dimension = &nc - 3;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let representation = cooking
            .try_cook(ColorAdjoint {}.to_symbolic([dimension.clone()]).as_view())
            .unwrap();
        let encoded = representation
            .as_fun_view()
            .unwrap()
            .iter()
            .next()
            .unwrap()
            .to_owned();
        assert_ne!(encoded, dimension);
        let prepared = |expression| PreparedNumerator {
            exact: NumeratorInterface::admit(expression, None).unwrap(),
            sample_source: NumeratorInterface::admit(Atom::one(), None).unwrap(),
            samples: Vec::new(),
        };
        let candidate = prepared(Atom::one() / encoded);
        let master = prepared(Atom::one());
        let scalars = scalar_names(&model());
        for concrete in [false, true] {
            let ratio = compare(
                &candidate,
                &master,
                ComparisonMode::UpToScalar,
                &GraphGroupingOptions {
                    check_canonical_numerator: true,
                    fully_numerical_substitution: concrete,
                    number_of_numerical_samples: 0,
                    ..Default::default()
                },
                &scalars,
            )
            .unwrap();
            if concrete {
                assert!(ratio.is_none());
            } else {
                assert_eq!(ratio, Some(Atom::one() / &dimension));
            }
        }
    }

    #[test]
    fn sampled_scalar_ratios_use_exact_common_denominators() {
        let model = model();
        let master = PreparedNumerator {
            exact: SymbolicTensor::infer(Atom::one()).unwrap(),
            sample_source: SymbolicTensor::infer(Atom::one()).unwrap(),
            samples: [test_atom("M-N"), Atom::num(3)]
                .map(|value| SymbolicTensor::infer(value).unwrap())
                .to_vec(),
        };
        let mut candidate = PreparedNumerator {
            exact: SymbolicTensor::infer(Atom::one()).unwrap(),
            sample_source: SymbolicTensor::infer(Atom::one()).unwrap(),
            samples: [test_atom("M^2-N^2"), test_atom("3*(M+N)")]
                .map(|value| SymbolicTensor::infer(value).unwrap())
                .to_vec(),
        };
        let options = GraphGroupingOptions::default();
        let ratio = compare(
            &candidate,
            &master,
            ComparisonMode::UpToScalar,
            &options,
            &scalar_names(&model),
        )
        .unwrap()
        .unwrap();
        assert!((ratio - test_atom("M+N")).together().is_zero());
        candidate.samples[1] = SymbolicTensor::infer(test_atom("3*(M+N+1)")).unwrap();
        assert!(
            compare(
                &candidate,
                &master,
                ComparisonMode::UpToScalar,
                &options,
                &scalar_names(&model)
            )
            .unwrap()
            .is_none()
        );
    }

    #[test]
    fn conjugated_model_couplings_are_scalar_but_conjugated_tensors_are_not() {
        let model = model();
        let scalars = scalar_names(&model);
        let m = test_atom("M");
        let n = test_atom("N");
        let common = symbol!("FeynKit::conjugated_ratio_tensor").call(1);
        for head in [Symbol::CONJ, INBUILTS.conj] {
            let master = &common * function!(head, &m);
            let candidate = &common * function!(head, &n) * 2;
            let ratio = compare_atoms(&candidate, &master, ComparisonMode::UpToScalar, &scalars)
                .unwrap()
                .expect("a ratio of conjugated model couplings is scalar");
            assert_eq!(ratio, 2 * function!(head, &n) / function!(head, &m));
            let tensor = test_atom("spenso::g(spenso::coad(8,a),spenso::coad(8,b))");
            assert!(!expression_is_scalar(&function!(head, tensor), &scalars));
        }
        assert!(!expression_is_scalar(
            &function!(INBUILTS.conj, 1, 2),
            &scalars
        ));
        assert!(!expression_is_scalar(
            &test_atom("foreign::conj(M)"),
            &scalars
        ));
    }

    #[test]
    fn canonical_only_comparison_honors_concrete_color_dimensions() {
        let model = Arc::new(model());
        let nc = Atom::var(CS.nc);
        for concrete in [false, true] {
            let grouped = group_diagrams(
                vec![
                    diagram(&model, "g0", "x", 25, 25),
                    diagram(&model, "g1", &format!("{}*x", nc.to_plain_string()), 25, 25),
                ]
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
                &model,
                &Grouping::UpToScalar(GraphGroupingOptions {
                    fully_numerical_substitution: concrete,
                    number_of_numerical_samples: 0,
                    ..exact_options()
                }),
                false,
                &Run::default(),
                false,
                DiagramLabels::default(),
            )
            .unwrap();
            assert_eq!(grouped.groups.len(), 1);
            assert_eq!(
                grouped.groups[0].members[1].ratio,
                if concrete { Atom::num(3) } else { nc.clone() }
            );
            assert_eq!(
                grouped.diagrams[0].numerator().unwrap().expression(),
                &test_atom("x")
            );
        }
    }

    #[test]
    fn canonical_scalar_rescaling_collects_only_numerical_sum_content() {
        let model = Arc::new(model());
        let master = diagram(&model, "master", "(M+N)^7*(1i/9*x-1i/9*y)", 25, 25);
        let source = master.numerator().unwrap().clone();
        let grouped = group_diagrams(
            vec![master, diagram(&model, "member", "(M+N)^7*(y-x)", 25, 25)]
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
            &model,
            &Grouping::UpToScalar(GraphGroupingOptions {
                number_of_numerical_samples: 0,
                ..exact_options()
            }),
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members[1].ratio, 9 * Atom::i());
        assert_eq!(grouped.diagrams[0].numerator().unwrap(), &source);
    }

    #[test]
    fn completed_samples_contract_closed_color_metric_coefficients() {
        let model = Arc::new(model());
        let diagram = diagram(&model, "color_loop", "1", 25, 25);
        let spectator = test_atom("(M+N)^7");
        for dimension in ["3", "spenso::Nc"] {
            let source = test_atom(&format!(
                "(M+N)^7
                 *spenso::g(spenso::dind(spenso::cof({dimension},FeynKit::DummyIndex(0))),spenso::cof({dimension},FeynKit::DummyIndex(1)))
                 *spenso::g(spenso::dind(spenso::cof({dimension},FeynKit::DummyIndex(1))),spenso::cof({dimension},FeynKit::DummyIndex(0)))"
            ));
            for sample in 0..3 {
                assert_eq!(
                    evaluate_tensor_sample(
                        &source,
                        &diagram,
                        &GraphGroupingOptions::default(),
                        sample
                    )
                    .unwrap(),
                    3 * &spectator
                );
            }
            assert!(source.contains(&spectator));
        }
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
        .unwrap()
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
            .unwrap()
            .is_none(),
            "a quotient retaining tensor heads must not be accepted as a scalar ratio"
        );
    }

    #[test]
    fn canonical_scalar_ratios_decode_dimensions_without_decoding_indices() {
        let model = Arc::new(model());
        let nc = Atom::var(CS.nc);
        let dimension = nc.clone().pow(2) - Atom::one();
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let representation = cooking
            .try_cook(ColorAdjoint {}.to_symbolic([dimension.clone()]).as_view())
            .unwrap();
        let encoded_dimension = representation
            .as_fun_view()
            .unwrap()
            .iter()
            .next()
            .unwrap()
            .to_owned();
        assert_ne!(encoded_dimension, dimension);

        for concrete in [false, true] {
            let grouped = group_diagrams(
                vec![
                    diagram(
                        &model,
                        "master",
                        &format!("{}*x", nc.to_canonical_string()),
                        25,
                        25,
                    ),
                    diagram(
                        &model,
                        "member",
                        &format!("{}*x", encoded_dimension.to_canonical_string()),
                        25,
                        25,
                    ),
                ]
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
                &model,
                &Grouping::UpToScalar(GraphGroupingOptions {
                    fully_numerical_substitution: concrete,
                    number_of_numerical_samples: 0,
                    ..exact_options()
                }),
                false,
                &Run::default(),
                false,
                DiagramLabels::default(),
            )
            .unwrap();
            assert_eq!(grouped.groups.len(), 1);
            assert_eq!(
                grouped.groups[0].members[1].ratio,
                if concrete {
                    Atom::num(8) / Atom::num(3)
                } else {
                    (&dimension / &nc).collect_factors()
                },
            );
        }

        let index = cooking
            .cook_function(
                function!(feynkit_graph::symbols::hedge_index(), 17, 2)
                    .as_fun_view()
                    .unwrap(),
            )
            .unwrap();
        assert!(index.get_symbol().is_some());
        let common = symbol!("FeynKit::dimension_ratio_tensor").call(1);
        assert!(
            compare_atoms(
                &(&index * &common),
                &common,
                ComparisonMode::UpToScalar,
                &scalar_names(&model),
            )
            .unwrap()
            .is_none(),
            "an encoded index must not be accepted as a scalar dimension coefficient",
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
            diagrams
                .clone()
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
            &model,
            &Grouping::Identical(exact_options()),
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
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
            diagrams.into_iter().map(DiagramPair::Shared).collect(),
            &model,
            &Grouping::Identical(exact_species),
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 4);
    }

    #[test]
    fn distinct_cut_metadata_is_merged_into_the_denominator_master() {
        let model = Arc::new(model());
        let diagrams = vec![cut_diagram(&model, false), cut_diagram(&model, true)];
        let ids = diagrams
            .iter()
            .map(|diagram| diagram.id().unwrap())
            .collect::<Vec<_>>();
        let source_cuts = diagrams
            .iter()
            .flat_map(|diagram| diagram.cuts().iter().cloned())
            .collect::<BTreeSet<_>>();
        let grouped = group_diagrams(
            diagrams.into_iter().map(DiagramPair::Shared).collect(),
            &model,
            &Grouping::Identical(exact_options()),
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
        )
        .unwrap();

        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members.len(), 2);
        assert_eq!(grouped.diagrams.len(), 1);
        // The staged set reads source identities from its ungrouped diagrams: the merged
        // cuts give the master a new identity.
        assert!(
            ids.iter()
                .all(|id| *id != grouped.diagrams[0].id().unwrap())
        );
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
            vec![master, member]
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
            &model,
            &Grouping::Identical(exact_options()),
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
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
            vec![master, member]
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
            &model,
            &Grouping::Identical(exact_options()),
            true,
            &Run::default(),
            false,
            DiagramLabels::default(),
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
            vec![identity.clone(), global]
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
            &model,
            &Grouping::Identical(exact_options()),
            true,
            &Run::default(),
            false,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members.len(), 2);

        let grouped = group_diagrams(
            vec![identity, partial]
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
            &model,
            &Grouping::Identical(exact_options()),
            true,
            &Run::default(),
            false,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 2);
    }

    #[test]
    fn removes_only_proven_zero_numerators_and_preserves_source_indices() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x-x", 25, 25),
                diagram(&model, "g1", "1", 25, 25),
            ]
            .into_iter()
            .map(DiagramPair::Shared)
            .collect(),
            &model,
            &Grouping::None,
            false,
            &Run::default(),
            true,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.zero_numerator_count, 1);
        assert_eq!(grouped.diagrams.len(), 1);
        assert_eq!(grouped.groups[0].members[0].diagram, 0);
        assert_eq!(grouped.groups[0].members[0].source_diagram, 1);
        assert_eq!(grouped.diagrams[0].name(), Some("g1"));
    }

    #[test]
    fn cancelling_weighted_groups_preserve_surviving_source_identity() {
        let model = Arc::new(model());
        let survivor = diagram(&model, "survivor", "(x+y)^7", 45, 25);
        let numerator = survivor.numerator().unwrap().clone();
        let grouped = group_diagrams(
            vec![
                diagram(&model, "half", "x+y", 25, 25).with_overall_factor(Atom::num((1, 2))),
                diagram(&model, "minus_quarter", "2*x+2*y", 25, 25)
                    .with_overall_factor(Atom::num((-1, 4))),
                survivor,
            ]
            .into_iter()
            .map(DiagramPair::Shared)
            .collect(),
            &model,
            &Grouping::UpToScalar(exact_options()),
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.zero_numerator_count, 0);
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.diagrams.len(), 1);
        assert_eq!(grouped.diagrams[0].name(), Some("survivor"));
        assert_eq!(grouped.diagrams[0].numerator().unwrap(), &numerator);
        assert_eq!(grouped.groups[0].master, 0);
        assert_eq!(grouped.groups[0].members[0].diagram, 0);
        assert_eq!(grouped.groups[0].members[0].source_diagram, 2);
    }

    #[test]
    fn grouped_scalar_weights_cancel_after_rational_rescaling() {
        let model = Arc::new(model());
        // The independently checked FullSM two-loop GL51--56 scalar class.
        // Its normalization vanishes only after combining rational weights.
        let ratio = test_atom(
            "(-5/3*UFO::cw^(-1)*UFO::sw+UFO::cw*UFO::sw^(-1))^(-1)*1/3*UFO::cw*UFO::sw^(-1)
            +(-5/3*UFO::cw^(-1)*UFO::sw+UFO::cw*UFO::sw^(-1))^(-1)*1/3*UFO::cw^(-1)*UFO::sw",
        );
        let opposite = test_atom(
            "(-5/3*UFO::cw^(-1)*UFO::sw+UFO::cw*UFO::sw^(-1))^(-1)*-1*UFO::cw*UFO::sw^(-1)
            +(-5/3*UFO::cw^(-1)*UFO::sw+UFO::cw*UFO::sw^(-1))^(-1)*1/3*UFO::cw^(-1)*UFO::sw",
        );
        let cancelling_classes = [
            vec![
                Atom::one(),
                ratio.clone(),
                ratio.clone(),
                ratio,
                opposite.clone(),
                opposite,
            ],
            // The charged-ghost/G0 tadpole normalization in the same model.
            vec![Atom::one(), test_atom("-(1/2)^(1/2)*2^(1/2)")],
        ];
        for (weights, scale) in cancelling_classes.iter().flat_map(|weights| {
            [Atom::one(), test_atom("UFO::ee^2/UFO::G^2")]
                .into_iter()
                .map(move |scale| (weights, scale))
        }) {
            for residual in [Atom::Zero, Atom::num((-1, 2))] {
                let diagrams = weights
                    .iter()
                    .chain([&residual])
                    .enumerate()
                    .map(|(index, weight)| {
                        diagram(&model, &format!("g{index}"), "(x+y)^7", 25, 25)
                            .with_overall_factor(weight * &scale)
                    })
                    .collect::<Vec<_>>();
                let grouped = group_diagrams(
                    diagrams.into_iter().map(DiagramPair::Shared).collect(),
                    &model,
                    &Grouping::UpToScalar(exact_options()),
                    false,
                    &Run::default(),
                    false,
                    DiagramLabels::default(),
                )
                .unwrap();
                if residual.is_zero() {
                    assert!(grouped.diagrams.is_empty());
                    assert!(grouped.groups.is_empty());
                } else {
                    assert_eq!(grouped.diagrams.len(), 1);
                    assert_eq!(
                        grouped.diagrams[0].numerator().unwrap().expression(),
                        &test_atom("(x+y)^7")
                    );
                    assert!(
                        evaluate_overall_factor(grouped.diagrams[0].overall_factor().as_view())
                            .check_equal_with(
                                &residual * &scale,
                                &symbolica::atom::ZeroTestSettings::default()
                                    .finite_field_samples(0)
                                    .numerical_samples(0)
                                    .exact_fallback(true),
                            )
                            .is_true()
                    );
                }
            }
        }
    }

    #[test]
    fn completed_samples_can_certify_a_factorized_scalar_zero() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram(&model, "zero", "(M+N)^2-M^2-2*M*N-N^2", 25, 25),
                diagram(&model, "survivor", "(M+N)^7", 45, 25),
            ]
            .into_iter()
            .map(DiagramPair::Shared)
            .collect(),
            &model,
            &Grouping::UpToScalar(GraphGroupingOptions {
                check_canonical_numerator: false,
                number_of_numerical_samples: 2,
                ..GraphGroupingOptions::default()
            }),
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.zero_numerator_count, 1);
        assert_eq!(grouped.diagrams.len(), 1);
        assert_eq!(grouped.diagrams[0].name(), Some("survivor"));
        assert_eq!(
            grouped.diagrams[0].numerator().unwrap().expression(),
            &test_atom("(M+N)^7")
        );
        assert_eq!(grouped.groups[0].members[0].source_diagram, 1);
    }

    #[test]
    fn unproved_distributive_zero_numerators_remain_factorized_and_retained() {
        let model = Arc::new(model());
        let input = diagram(&model, "unproved-zero", "(x+y)^2-x^2-2*x*y-y^2", 25, 25);
        let numerator = input.numerator().unwrap().clone();
        assert!(!numerator.is_zero());
        let grouped = group_diagrams(
            vec![input].into_iter().map(DiagramPair::Shared).collect(),
            &model,
            &Grouping::None,
            false,
            &Run::default(),
            true,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.zero_numerator_count, 0);
        assert_eq!(grouped.diagrams.len(), 1);
        assert_eq!(grouped.diagrams[0].numerator().unwrap(), &numerator);
        assert_eq!(grouped.groups[0].members[0].source_diagram, 0);
    }

    #[test]
    fn zero_projectors_are_removed_before_grouping() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![diagram_with_projector(&model, "zero", "1", "0", 25, 25)]
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
            &model,
            &Grouping::None,
            false,
            &Run::default(),
            true,
            DiagramLabels::default(),
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
            vec![diagram(&model, "coupling-zero", "UFO::GC-1", 25, 25)]
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
            &model,
            &Grouping::None,
            false,
            &Run::default(),
            true,
            DiagramLabels::default(),
        )
        .unwrap();

        assert_eq!(grouped.zero_numerator_count, 1);
        assert!(grouped.diagrams.is_empty());
    }

    #[test]
    fn zero_detection_reduces_color_and_preserves_factorized_survivors() {
        let model = Arc::new(model());
        let [a, b, c, d, e] =
            [1, 2, 3, 4, 5].map(|index| ColorAdjoint {}.new_rep(8).to_symbolic([Atom::num(index)]));
        let bubble = idenso::color_f!(&a, &b, &c) * idenso::color_f!(&a, &b, &d);
        let zero = &bubble * idenso::color_f!(&c, &d, &e);
        assert!(!zero.is_zero());
        let spectator = test_atom("(x+y)^30");
        let inputs = vec![
            diagram_with_atoms(&model, "color-zero", zero * &spectator, Atom::one(), 25, 25),
            diagram_with_atoms(&model, "nonzero", bubble * &spectator, Atom::one(), 25, 25),
            diagram(&model, "unproved", "(x+y)^2-x^2-2*x*y-y^2", 25, 25),
        ];
        let expected = inputs[1..]
            .iter()
            .map(|diagram| diagram.numerator().unwrap().clone())
            .collect::<Vec<_>>();
        let unfiltered = group_diagrams(
            inputs
                .clone()
                .into_iter()
                .map(DiagramPair::Shared)
                .collect(),
            &model,
            &Grouping::None,
            false,
            &Run::default(),
            false,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(unfiltered.zero_numerator_count, 0);
        assert_eq!(unfiltered.diagrams.len(), inputs.len());
        for (original, retained) in inputs.iter().zip(&unfiltered.diagrams) {
            assert_eq!(original.to_json().unwrap(), retained.to_json().unwrap());
        }
        let grouped = group_diagrams(
            inputs.into_iter().map(DiagramPair::Shared).collect(),
            &model,
            &Grouping::None,
            false,
            &Run::default(),
            true,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.zero_numerator_count, 1);
        assert_eq!(
            grouped
                .diagrams
                .iter()
                .map(|diagram| diagram.numerator().unwrap().clone())
                .collect::<Vec<_>>(),
            expected
        );
        assert_eq!(
            grouped
                .groups
                .iter()
                .map(|group| group.members[0].source_diagram)
                .collect::<Vec<_>>(),
            vec![1, 2]
        );
    }

    #[test]
    fn zero_detection_combines_reduced_color_coefficients_and_tensor_powers() {
        let model = Arc::new(model());
        let [a, b, c, d] =
            [1, 2, 3, 4].map(|index| ColorAdjoint {}.new_rep(8).to_symbolic([Atom::num(index)]));
        let bubble = idenso::color_f!(&a, &b, &c) * idenso::color_f!(&a, &b, &d);
        let metric = function!(spenso::network::library::symbolic::ETS.metric, &c, &d);
        let lorentz = Minkowski {}.new_rep(4);
        let lorentz_metric = function!(
            spenso::network::library::symbolic::ETS.metric,
            lorentz.to_symbolic([Atom::num(1)]),
            lorentz.to_symbolic([Atom::num(2)])
        );
        let spectator = test_atom("(x+y)^30") * &lorentz_metric;
        // These zeros need actual color identities, beyond odd graph symmetries.
        // Different Lorentz coefficients of one color structure must combine
        // after reduction, while repeated color tensors must not be skipped.
        let mixed = &bubble * &spectator - Atom::num(3) * &metric * &spectator;
        // The coefficient sum and its individual terms receive different
        // aliases. Recollection must still expose their cancellation when the
        // definitions are resolved after the color rewrite.
        let other = test_atom("(u+v)^30") * lorentz_metric;
        let aliased =
            &bubble / Atom::num(3) * (&spectator + &other) - &metric * &spectator - &metric * other;
        let power = idenso::color_f!(&a, &b, &c).pow(Atom::num(2)) - Atom::num(24);
        let nonzero = &bubble * &spectator - Atom::num(2) * &metric * &spectator;
        let inputs = vec![
            diagram_with_atoms(&model, "mixed-zero", mixed, Atom::one(), 25, 25),
            diagram_with_atoms(&model, "power-zero", power, Atom::one(), 25, 25),
            diagram_with_atoms(&model, "alias-zero", aliased, Atom::one(), 25, 25),
            diagram_with_atoms(&model, "nonzero", nonzero.clone(), Atom::one(), 25, 25),
        ];
        let grouped = group_diagrams(
            inputs.into_iter().map(DiagramPair::Shared).collect(),
            &model,
            &Grouping::None,
            false,
            &Run::default(),
            true,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(grouped.zero_numerator_count, 3);
        assert_eq!(grouped.diagrams.len(), 1);
        assert_eq!(
            grouped.diagrams[0].numerator().unwrap().expression(),
            &nonzero
        );
        assert_eq!(grouped.groups[0].members[0].source_diagram, 3);
    }

    #[test]
    fn zero_detection_does_not_depend_on_index_labels() {
        let model = Arc::new(model());
        // Index payloads are cooked to symbols in the order a process first
        // meets them; for generated half-edge labels that order changes from
        // process to process. Relabelling the dummies of a color zero must not
        // change its proof, also beside an indexed momentum spectator. Each
        // copy uses fresh labels, cooked in its own term order.
        let numerator = |labels: [i64; 6]| {
            let index = |label: i64| function!(symbol!("feynkit_grouping_test::label"), label);
            let [a, b, c, d, e] = [0, 1, 2, 3, 4].map(|position| {
                ColorAdjoint {}
                    .new_rep(8)
                    .to_symbolic([index(labels[position])])
            });
            let lorentz = Minkowski {}.new_rep(4).to_symbolic([index(labels[5])]);
            let momentum = |edge: i64| {
                FunctionBuilder::new(kinematics_symbols::momentum())
                    .add_arg(edge)
                    .add_arg(&lorentz)
                    .finish()
            };
            idenso::color_f!(&a, &b, &c)
                * idenso::color_f!(&a, &b, &d)
                * idenso::color_f!(&c, &d, &e)
                * momentum(0)
                * momentum(2)
        };
        let inputs = (0..6)
            .flat_map(|shift| [(shift, false), (shift, true)])
            .map(|(shift, reversed)| {
                // Rotations and reflections of one label order.
                let offset = 1000 + 10 * (2 * shift + i64::from(reversed));
                let labels = [0, 1, 2, 3, 4, 5].map(|position| {
                    let position = (position + shift) % 6;
                    offset + if reversed { 5 - position } else { position }
                });
                diagram_with_atoms(
                    &model,
                    &format!("relabelled-{shift}-{reversed}"),
                    numerator(labels),
                    Atom::one(),
                    25,
                    25,
                )
            })
            .collect::<Vec<_>>();
        let grouped = group_diagrams(
            inputs.into_iter().map(DiagramPair::Shared).collect(),
            &model,
            &Grouping::None,
            false,
            &Run::default(),
            true,
            DiagramLabels::default(),
        )
        .unwrap();
        assert_eq!(
            grouped
                .diagrams
                .iter()
                .map(|diagram| diagram.name().unwrap())
                .collect::<Vec<_>>(),
            Vec::<&str>::new()
        );
        assert_eq!(grouped.zero_numerator_count, 12);
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn parallel_zero_detection_preserves_input_order_and_progress() {
        let model = Arc::new(model());
        let inputs = (0..32)
            .map(|index| {
                diagram(
                    &model,
                    &format!("g{index}"),
                    if index % 3 == 0 { "0" } else { "(x+y)^30" },
                    25,
                    25,
                )
            })
            .collect::<Vec<_>>();
        let expected_sources = (0..32).filter(|index| index % 3 != 0).collect::<Vec<_>>();
        for threads in [1, 4] {
            let updates = Arc::new(Mutex::new(Vec::new()));
            let captured = updates.clone();
            let options = Run::default().progress(move |snapshot| {
                captured.lock().unwrap().push(snapshot);
                crate::GenerationControl::Continue
            });
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let grouped = pool
                .install(|| {
                    group_diagrams(
                        inputs
                            .clone()
                            .into_iter()
                            .map(DiagramPair::Shared)
                            .collect(),
                        &model,
                        &Grouping::None,
                        false,
                        &options,
                        true,
                        DiagramLabels::default(),
                    )
                })
                .unwrap();
            assert_eq!(grouped.zero_numerator_count, 11);
            assert_eq!(
                grouped
                    .groups
                    .iter()
                    .map(|group| group.members[0].source_diagram)
                    .collect::<Vec<_>>(),
                expected_sources
            );
            for (diagram, &source) in grouped.diagrams.iter().zip(&expected_sources) {
                assert_eq!(diagram.name(), inputs[source].name());
                assert_eq!(
                    diagram.numerator().unwrap(),
                    inputs[source].numerator().unwrap()
                );
            }
            assert_eq!(
                updates
                    .lock()
                    .unwrap()
                    .iter()
                    .map(|progress| (
                        progress.stage,
                        usize::try_from(&progress.completed).unwrap(),
                        progress.total.as_ref().map(|n| usize::try_from(n).unwrap())
                    ))
                    .collect::<Vec<_>>(),
                (0..=32)
                    .map(|count| (crate::Stage::ZeroColor, count, Some(32)))
                    .collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn neutral_vector_forward_numerator_groups_chiral_vertices() {
        let mut model =
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
        model.simplify_zero_values().unwrap();
        let model = Arc::new(model);
        let process = Process::new([23_i64], [1_i64, -1]);
        let run = Run::default().threads(1);
        let generate = |grouping: Grouping| {
            process
                .clone()
                .cross_section(crate::Loops::forward(1..=1))
                .with_topology(crate::TopologyRules::new().with_max_vertices(2))
                .generate(Arc::clone(&model), &run)
                .unwrap()
                .build(
                    &crate::NumeratorOptions::default().with_naming(crate::Naming::prefix("GL")),
                    &run,
                )
                .unwrap()
                .group(&grouping, &run)
                .unwrap()
        };
        let raw = generate(Grouping::None);
        let grouped = generate(Grouping::UpToScalar(Default::default()));
        assert_eq!(raw.diagrams().len(), 1);
        assert_eq!(grouped.diagrams().len(), 1);
        assert_eq!(
            grouped.diagrams()[0].numerator().unwrap(),
            raw.diagrams()[0].numerator().unwrap()
        );
        assert_eq!(
            grouped.diagrams()[0].projector(),
            raw.diagrams()[0].projector()
        );
        assert_eq!(grouped.diagrams()[0].cuts(), raw.diagrams()[0].cuts());
    }
    #[test]
    #[ignore = "explicit WP1 original-pair certificate and production fallback matrix"]
    fn vacuum_group28_original_pair_identity() {
        let mut model = Model::standard_model();
        let card = ParameterCard::from_path(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../assets/models/json/sm/restrict_default.json"
        ))
        .unwrap();
        model.apply_parameter_card(&card).unwrap();
        model.simplify_zero_values().unwrap();
        let model = Arc::new(model);
        let dots = [
            include_str!("../tests/fixtures/vacuum_b2/GL028.dot"),
            include_str!("../tests/fixtures/vacuum_b2/GL034.dot"),
        ];
        // Frozen provenance is kept byte-for-byte. #209 moved only Q's q
        // presentation tag to notation; bind that known spelling to its current
        // constructor before parsing. No indices, coefficients or tensor
        // attributes are inferred or changed by this fixture-only adaptation.
        let legacy_q = "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q";
        let native_q = Atom::var(kinematics_symbols::momentum()).to_canonical_string();
        assert_eq!(
            native_q,
            "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::Q"
        );
        let legacy_edge = "gammalooprs::{spenso::index,spenso::index-label:e}::edge";
        let native_edge = Atom::var(graph_api::symbols::edge_index()).to_canonical_string();
        assert_eq!(native_edge, "gammalooprs::{spenso::index}::edge");
        let native_dots = dots.map(|dot| {
            assert!(dot.contains(legacy_q));
            assert!(dot.contains(legacy_edge));
            dot.replace(legacy_q, &native_q)
                .replace(legacy_edge, &native_edge)
        });
        let dots = native_dots.each_ref().map(String::as_str);
        // This finite association is pinned for these two frozen inputs only.
        // The production importer still rejects their historical Model identity.
        let historical = "e845ddde3d54cde80542b9c73ca6392b518dcabe8c6b5a551343e5f933ac7f98";
        let associated = "a9f6ae4fafba3163992d30b6fa9152860a351f0a97dbe1aeb74f499dc2d7df69";
        let current = model.fingerprint().to_string();
        assert_eq!(current, associated, "B2 Model association needs review");
        assert!(matches!(
            FeynmanDiagram::from_dot(Arc::clone(&model), dots[0]),
            Err(DiagramError::DotModelFingerprintMismatch { serialized, actual })
                if serialized == historical && actual.to_string() == current
        ));
        let historical_tag = format!("model_fingerprint=\"{historical}\"");
        let associated_tag = format!("model_fingerprint=\"{associated}\"");
        let associate_frozen = |dot: &str, target: &str| {
            (target == associated
                && dots.contains(&dot)
                && dot.matches(&historical_tag).count() == 1)
                .then(|| dot.replacen(&historical_tag, &associated_tag, 1))
        };
        assert!(associate_frozen(dots[0], historical).is_none());
        let foreign = dots[0].replacen("interaction_id=136", "interaction_id=135", 1);
        assert_ne!(foreign, dots[0]);
        assert!(associate_frozen(&foreign, &current).is_none());
        let missing = dots[0].replacen(&historical_tag, "model_fingerprint=\"foreign\"", 1);
        assert!(associate_frozen(&missing, &current).is_none());
        let duplicate = format!("{}\n{historical_tag}", dots[0]);
        assert!(associate_frozen(&duplicate, &current).is_none());
        let associated_dots = dots.map(|dot| {
            associate_frozen(dot, &current).expect("verified frozen B2 input association")
        });
        let dots = associated_dots.each_ref().map(String::as_str);
        let diagrams = dots.map(|dot| FeynmanDiagram::from_dot(Arc::clone(&model), dot).unwrap());
        let originals = diagrams.each_ref().map(|diagram| {
            // These frozen sources have unit projectors and prefactors; assert
            // that precondition instead of invoking any extra tensor algebra.
            assert!(diagram.numerator_prefactor().expression().is_one());
            assert!(diagram.projector().expression().is_one());
            let complete = diagram.numerator().unwrap();
            let complete = complete
                .with_checked_expression(model.expand_couplings(complete.expression()))
                .unwrap();
            super::normalize_momentum_routing(diagram, &complete).unwrap()
        });
        // Only color heads are selected: the propagator sums are immutable
        // Lorentz spectators. Reconstruct the source product syntactically.
        let library: GroupingTensorLibrary = hep_lib::hep_lib_atom();
        let functions = Panic::new_lib();
        let parsing = TensorParseSettings {
            shorthand_parsing: ShorthandParsing::expand_schoonschip_only(),
            ..Default::default()
        };
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let sectors = originals.each_ref().map(|source| {
            let mut sectors = source
                .coefficient_list(TensorCollectFilter::Reps([
                    ColorAdjoint {}.into(),
                    ColorFundamental {}.into(),
                    ColorSextet {}.into(),
                ]))
                .unwrap();
            assert_eq!(sectors.len(), 1);
            let (color, lorentz) = sectors.pop().unwrap();
            assert_eq!(
                (color.expression() * lorentz.expression()).collect_factors(),
                source.expression().collect_factors()
            );
            let lowered: Canonicalized<GroupingTensorNetwork> = color
                .into_network(&library, &functions, &parsing, Some(&cooking))
                .unwrap();
            let (mut network, layout) = lowered.into_parts();
            assert_eq!(layout.order(), 0);
            network
                .execute::<Sequential, MinResultRank, _, _, _>(&library, &functions)
                .unwrap();
            let color: Atom = network.result_scalar().unwrap().into();
            (color.together(), lorentz)
        });
        eprintln!(
            "isolated finite SU(3) color factors: {}, {}",
            sectors[0].0, sectors[1].0
        );

        assert_eq!(sectors[0].0, Atom::num(7) / Atom::num(3));
        assert_eq!(sectors[1].0, -Atom::num(2) / Atom::num(3));

        // Verify C^2=-I and C gamma^T C=gamma using exact finite matrix
        // components, with every Lorentz component checked (no momentum point).
        let gamma_key = AGS.gamma_strct::<GroupingIndex>(4);
        let gamma_layout = TensorDataLayout::from_canonicalized(&gamma_key).unwrap();
        let gamma = hep_lib::gamma_data_weyl(gamma_key, 1_i32, 0_i32);
        let component = |row, column, mu| {
            let indices = gamma_layout
                .logical_expanded_to_storage_expanded(&[row, column, mu])
                .unwrap();
            let value = gamma.canonical().smart_get(&indices).unwrap();
            Atom::num(value.re) + Atom::i() * Atom::num(value.im)
        };
        let charge = |row, column| {
            CHARGE_CONJUGATION_WEYL_COMPONENTS
                .iter()
                .find_map(|(indices, sign)| (*indices == [row, column]).then_some(Atom::num(*sign)))
                .unwrap_or(Atom::Zero)
        };
        for row in 0..4 {
            for column in 0..4 {
                assert_eq!(charge(row, column), -charge(column, row));
                let square = Atom::add_many((0..4).map(|a| charge(row, a) * charge(a, column)));
                assert_eq!(
                    square,
                    if row == column {
                        -Atom::one()
                    } else {
                        Atom::Zero
                    }
                );
                for mu in 0..4 {
                    let sandwich = Atom::add_many((0..4).flat_map(|a| {
                        (0..4)
                            .map(move |b| charge(row, a) * component(b, a, mu) * charge(b, column))
                    }));
                    assert_eq!(sandwich, component(row, column, mu));
                }
            }
        }
        // Vertices 1,2,3 and edges 3,4,6 form one closed top loop.
        // Applying the same similarity C^-1 A C to every local matrix
        // preserves its trace and never opens a propagator sum.
        let loop_vertices = BTreeSet::from([VertexId(1), VertexId(2), VertexId(3)]);
        let mut spinor_slots = BTreeSet::new();
        let mut loop_vertex_gammas = 0;
        for (id, vertex) in diagrams[0].vertices() {
            if !loop_vertices.contains(&id) {
                continue;
            }
            vertex.numerator.expression().visitor(&mut |part| {
                if let AtomView::Fun(matrix) = part
                    && matrix.get_symbol() == AGS.gamma
                {
                    loop_vertex_gammas += 1;
                    spinor_slots.extend(matrix.iter().take(2).map(|slot| slot.to_owned()));
                }
                true
            });
        }
        assert_eq!(loop_vertex_gammas, 3);
        assert_eq!(spinor_slots.len(), 6);
        let loop_edges = diagrams[0]
            .edges()
            .filter(|(_, _, edge)| model.particle_by_id(edge.particle).unwrap().is_fermion())
            .filter(|(_, endpoints, _)| {
                endpoints.source.is_some_and(|v| loop_vertices.contains(&v))
                    || endpoints.target.is_some_and(|v| loop_vertices.contains(&v))
            })
            .map(|(id, endpoints, _)| {
                assert!(endpoints.source.is_some_and(|v| loop_vertices.contains(&v)));
                assert!(endpoints.target.is_some_and(|v| loop_vertices.contains(&v)));
                id
            })
            .collect::<BTreeSet<_>>();
        assert_eq!(
            loop_edges,
            BTreeSet::from([EdgeId(3), EdgeId(4), EdgeId(6)])
        );
        // Audit the similarity's entire matrix domain: only four-dimensional
        // gammas and identities occur, including all three mass alternatives.
        let mut loop_identities = 0;
        let mut gamma_slot_uses = BTreeMap::<Atom, usize>::new();
        let mut identity_slot_uses = BTreeMap::<Atom, usize>::new();
        sectors[0].1.expression().visitor(&mut |part| {
            if let AtomView::Fun(matrix) = part {
                let args = matrix.iter().collect::<Vec<_>>();
                if args
                    .iter()
                    .any(|arg| spinor_slots.contains(&arg.to_owned()))
                {
                    assert!(spinor_slots.contains(&args[0].to_owned()));
                    assert!(spinor_slots.contains(&args[1].to_owned()));
                    if matrix.get_symbol() == AGS.gamma {
                        for slot in &args[..2] {
                            *gamma_slot_uses.entry(slot.to_owned()).or_default() += 1;
                        }
                        assert_eq!(args.len(), 3);
                        let AtomView::Fun(lorentz_slot) = args[2] else {
                            panic!("gamma Lorentz slot")
                        };
                        assert_eq!(lorentz_slot.iter().next().unwrap(), Atom::num(4).as_view());
                    } else {
                        assert_eq!(matrix.get_symbol(), ETS.metric);
                        assert_eq!(args.len(), 2);
                        loop_identities += 1;
                        for slot in &args[..2] {
                            *identity_slot_uses.entry(slot.to_owned()).or_default() += 1;
                        }
                    }
                }
            }
            true
        });
        assert_eq!(loop_identities, 3);
        for slot in &spinor_slots {
            assert_eq!(gamma_slot_uses[slot], 2);
            assert_eq!(identity_slot_uses[slot], 1);
            let AtomView::Fun(slot) = slot.as_view() else {
                panic!("bispinor slot")
            };
            assert_eq!(slot.iter().next().unwrap(), Atom::num(4).as_view());
        }
        let mut transformed_gammas = 0;
        let transformed = sectors[0].1.expression().replace_map(|part, _, out| {
            let AtomView::Fun(matrix) = part else {
                return;
            };
            if matrix.get_symbol() != AGS.gamma {
                return;
            }
            let args = matrix.iter().collect::<Vec<_>>();
            if !spinor_slots.contains(&args[0].to_owned()) {
                return;
            }
            assert!(spinor_slots.contains(&args[1].to_owned()));
            transformed_gammas += 1;
            **out = -FunctionBuilder::new(AGS.gamma)
                .add_arg(args[1])
                .add_arg(args[0])
                .add_arg(args[2])
                .finish();
        });
        assert_eq!(transformed_gammas, 6);
        let normalize = |expression: Atom| {
            let expression = expression
                .replace_map_bottom_up(|part, _, out| {
                    if matches!(part, AtomView::Add(_)) {
                        **out = part.collect_numeric_content();
                    }
                })
                .collect_factors();
            // Pure graph/dummy relabelling: no tensor finisher, contraction,
            // Dirac identities, or algebra prerequisites are enabled here.
            expression
                .canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap()
        };
        let transformed = normalize(transformed);
        let master = normalize(-sectors[1].1.expression());
        assert_eq!(
            transformed, master,
            "closed-loop similarity must give L28=-L34 without distributing source factors"
        );
        let ratio = -(&sectors[0].0 / &sectors[1].0).together();
        assert_eq!(ratio, Atom::num(7) / Atom::num(2));
        eprintln!("factor-preserving original-pair identity: N28=(7/2)*N34");
        let sample_sources = sectors.map(|(color, lorentz)| {
            lorentz
                .with_checked_expression(lorentz.expression() * color)
                .unwrap()
        });
        let exact_only = GraphGroupingOptions {
            check_canonical_numerator: true,
            number_of_numerical_samples: 0,
            ..Default::default()
        };
        let source_indices = [34, 28];
        let exact_inputs = dots
            .iter()
            .rev()
            .map(|dot| {
                let diagram = FeynmanDiagram::from_dot(Arc::clone(&model), dot).unwrap();
                DiagramPair::Shared(diagram)
            })
            .collect();
        let exact_grouped = group_diagrams(
            exact_inputs,
            &model,
            &Grouping::UpToScalar(exact_only),
            false,
            &Run::default().threads(1),
            false,
            DiagramLabels {
                source_indices: Some(&source_indices),
                ..Default::default()
            },
        )
        .unwrap();
        // The actual frozen pair is not merged by the production exact path.
        // Pin both singleton owners, their unit ratios and their scalar weights;
        // a dropped source is not a successful zero-sample control.
        assert_eq!(exact_grouped.zero_numerator_count, 0);
        assert_eq!(exact_grouped.diagrams.len(), 2);
        assert_eq!(exact_grouped.groups.len(), 2);
        for (index, group) in exact_grouped.groups.iter().enumerate() {
            assert_eq!(group.master, index);
            assert_eq!(group.members.len(), 1);
            let member = &group.members[0];
            assert_eq!(member.source_diagram, source_indices[index]);
            assert_eq!(member.diagram, index);
            assert_eq!(member.ratio, Atom::one());
            let original = &diagrams[1 - index];
            assert_eq!(&member.overall_factor, original.overall_factor());
            assert_eq!(exact_grouped.diagrams[index].name(), original.name());
            assert_eq!(
                exact_grouped.diagrams[index].numerator().unwrap(),
                original.numerator().unwrap()
            );
            assert_eq!(
                evaluate_overall_factor(exact_grouped.diagrams[index].overall_factor().as_view()),
                evaluate_overall_factor(original.overall_factor().as_view())
            );
        }
        assert_eq!(
            exact_grouped
                .groups
                .iter()
                .flat_map(|group| group.members.iter().map(|member| member.source_diagram))
                .collect::<BTreeSet<_>>(),
            BTreeSet::from(source_indices)
        );
        eprintln!(
            "production samples=0 retains both original singleton owners; samples=5 below must merge through fallback"
        );
        for seed in [3, 7, 11] {
            for canonical in [false, true] {
                let options = GraphGroupingOptions {
                    numerical_sample_seed: seed,
                    check_canonical_numerator: canonical,
                    ..Default::default()
                };
                // Feed original local factors through the production owner,
                // master first, matching the GL034 -> GL028 direction.
                let inputs = dots
                    .iter()
                    .rev()
                    .map(|dot| {
                        let diagram = FeynmanDiagram::from_dot(Arc::clone(&model), dot).unwrap();
                        DiagramPair::Shared(diagram)
                    })
                    .collect();
                let grouping = Grouping::UpToScalar(options);
                let grouped = group_diagrams(
                    inputs,
                    &model,
                    &grouping,
                    false,
                    &Run::default().threads(1),
                    false,
                    DiagramLabels {
                        source_indices: Some(&source_indices),
                        ..Default::default()
                    },
                )
                .unwrap();
                assert_eq!(
                    grouped.groups.len(),
                    1,
                    "seed={seed}, canonical={canonical}"
                );
                assert_eq!(grouped.groups[0].members.len(), 2);
                assert_eq!(grouped.groups[0].master, 0);
                assert_eq!(grouped.zero_numerator_count, 0);
                assert_eq!(grouped.diagrams.len(), 1);
                assert_eq!(
                    grouped.groups[0]
                        .members
                        .iter()
                        .map(|member| member.source_diagram)
                        .collect::<Vec<_>>(),
                    source_indices
                );
                assert_eq!(grouped.groups[0].members[0].ratio, Atom::one());
                for member in &grouped.groups[0].members {
                    assert_eq!(member.diagram, 0);
                    assert_eq!(
                        &member.overall_factor,
                        diagrams[match member.source_diagram {
                            28 => 0,
                            34 => 1,
                            _ => panic!("unknown source identity"),
                        }]
                        .overall_factor()
                    );
                }
                let expected_weight =
                    evaluate_overall_factor(diagrams[1].overall_factor().as_view())
                        + &ratio * evaluate_overall_factor(diagrams[0].overall_factor().as_view());
                assert_eq!(expected_weight, Atom::num(3) / Atom::num(4));
                assert_eq!(
                    evaluate_overall_factor(grouped.diagrams[0].overall_factor().as_view())
                        .together(),
                    expected_weight
                );

                let stored = &grouped.groups[0].members[1].ratio;
                let Grouping::UpToScalar(options) = &grouping else {
                    unreachable!()
                };
                assert_eq!(
                    stored.together(),
                    ratio,
                    "seed={seed}, canonical={canonical}"
                );
                // Check the accepted stored ratio against both original sources at unused points.
                for sample in [7, 13] {
                    let candidate = super::evaluate_tensor_sample(
                        &sample_sources[0],
                        &diagrams[0],
                        options,
                        sample,
                    )
                    .unwrap();
                    let master = super::evaluate_tensor_sample(
                        &sample_sources[1],
                        &diagrams[1],
                        options,
                        sample,
                    )
                    .unwrap();
                    assert!(
                        PreparedNumerator::finite_samples_equal(
                            candidate.expression(),
                            &(stored.together() * master.expression())
                        ),
                        "seed={seed}, sample={sample}, canonical={canonical}"
                    );
                }
            }
        }
    }
}
