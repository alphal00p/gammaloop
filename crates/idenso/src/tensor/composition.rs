//! Port planning and composition of admitted symbolic tensors.
//!
//! See the [materialization](SymbolicTensor::materialized) and
//! [matrix composition](SymbolicTensor::compose) examples.

use super::inference::InterfaceInference;
use super::{SymbolicTensor, TensorProofs};
use crate::tensor::inference::TensorInferenceError;
use crate::{CookSettings, CookingError};
use ahash::AHashMap;
use linnet::half_edge::subgraph::ModifySubSet;
use linnet::half_edge::subgraph::subset::SubSet;
use spenso::network::library::symbolic::ETS;
use spenso::network::parsing::{
    AtomStructureExt, ChainNestingError, ParseState, StrictTensorFilter,
};
use spenso::network::tags::{SPENSO_TAG, SpensoTags};
use spenso::shadowing;
use spenso::shadowing::{ANTISYM, CYCLIC, SYM};
use spenso::structure::abstract_index::{AIND_SYMBOLS, AbstractIndex};
use spenso::structure::dimension::Dimension;
use spenso::structure::partial::{
    PartialIndex, PartialSlot, PartialStructure, PartialStructureExt,
};
use spenso::structure::representation::{LibraryRep, LibrarySlot, RepName, Representation};
use spenso::structure::slot::{IsAbstractSlot, Slot, SlotMatch, SlotMatcher};
use spenso::structure::{OrderedStructure, SlotIndex, TensorStructure};
use std::cell::OnceCell;
use std::cmp::Ordering as CmpOrdering;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::ops::Mul;
use std::sync::{Arc, OnceLock};
use std::{iter, mem, slice};
use symbolica::atom::representation::FunView;
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol};
use symbolica::coefficient::CoefficientView;
use symbolica::domains::rational::Rational;
use thiserror::Error;

mod chain;

#[cfg(test)]
thread_local! {
    pub(crate) static DUMMY_RESERVATION_WALKS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

impl SymbolicTensor<PartialStructure> {
    /// Build raw storage and finish product spelling. This does not admit the
    /// declared interface; external parts enter through `from_declared_parts`.
    /// Move existing carriers or use their typed operations to retain metadata.
    /// See the [admission example](SymbolicTensor).
    pub fn new(atom: Atom, interface: PartialStructure) -> Self {
        let atom = match Self::normalize_products(atom.as_view()) {
            AtomOrView::View(_) => atom,
            normalized => normalized.into_owned(),
        };
        Self::from_normalized_parts(atom, interface)
    }

    /// Attach raw parts without finishing or admitting them. This starts with
    /// no index history and cannot replace a metadata-preserving transfer of an
    /// existing tensor. Typed algebra and storage bindings retain the carrier.
    /// See the [admission example](SymbolicTensor).
    pub fn from_normalized_parts(atom: Atom, interface: PartialStructure) -> Self {
        let is_metric =
            matches!(atom.as_view(), AtomView::Fun(fun) if fun.get_symbol() == ETS.metric);
        let is_composite = !matches!(atom.as_view(), AtomView::Fun(_));
        Self {
            proofs: Default::default(),
            expression: atom,
            structure: interface.canonicalize_open_ports(),
            is_metric,
            is_composite,
        }
    }

    pub fn rank(&self) -> usize {
        self.structure.canonical().order()
    }

    pub fn is_scalar(&self) -> bool {
        self.rank() == 0
    }

    /// Exact normalized tensor-value equality, including public logical port order.
    /// Derived implementation flags, provenance, carried index history and proofs
    /// are not mathematical data. Value equality cannot certify identical future
    /// collision checks or interchangeable cached facts.
    pub fn same_value(&self, other: &Self) -> bool {
        self.expression == other.expression
            && self.structure.logical_slots() == other.structure.logical_slots()
    }

    /// Hash the same mathematical data used by `same_value`.
    pub fn hash_value<H: Hasher>(&self, state: &mut H) {
        self.expression.hash(state);
        self.structure.logical_slots().hash(state);
    }

    /// Whether the normalized symbolic payload is exactly zero.
    /// This does not expand the expression or erase the interface of a typed zero.
    pub fn is_zero(&self) -> bool {
        self.expression.is_zero()
    }

    /// Brackets are associative ordered products, not index scopes. Splice only
    /// bracket children; sums and chain/trace binders remain intact.
    fn bracket_product(left: AtomView<'_>, right: AtomView<'_>) -> Atom {
        let product = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(left)
            .add_arg(right)
            .finish();
        Self::normalize_products(product.as_view()).into_owned()
    }

    /// Ordinary multiplication is safe once tensor slots have identities.
    /// An unresolved occurrence must survive separately and in positional order:
    /// in particular, `bracket(p(rep), p(rep))` must not become `p(rep)^2`.
    pub(super) fn has_unresolved_ports(value: AtomView<'_>) -> bool {
        if let Ok(slot) = Slot::<LibraryRep, AbstractIndex>::try_from(value) {
            return matches!(slot.aind(), AbstractIndex::Open { .. });
        }
        if Representation::<LibraryRep>::try_from(value).is_ok() {
            return true;
        }
        match value {
            AtomView::Add(sum) => sum.iter().any(Self::has_unresolved_ports),
            AtomView::Mul(product) => product.iter().any(Self::has_unresolved_ports),
            AtomView::Pow(power) => Self::has_unresolved_ports(power.get_base_exp().0),
            // Dot owns the slots of its operands; they are not public occurrences.
            AtomView::Fun(fun)
                if fun.get_symbol() == SPENSO_TAG.dot || fun.get_symbol().is_scalar() =>
            {
                false
            }
            AtomView::Fun(fun)
                if is_composite_head(fun.get_symbol(), &SPENSO_TAG)
                    || fun.get_symbol() == AIND_SYMBOLS.aind =>
            {
                fun.iter()
                    .skip(usize::from(fun.get_symbol() == SPENSO_TAG.trace))
                    .any(Self::has_unresolved_ports)
            }
            AtomView::Fun(fun) if value.is_tensorial(StrictTensorFilter::Tagged) => fun
                .iter()
                .filter(|&arg| {
                    direct_structural_port(arg)
                        || matches!(arg, AtomView::Fun(f) if f.get_symbol() == AIND_SYMBOLS.aind)
                })
                .any(Self::has_unresolved_ports),
            _ => false,
        }
    }

    fn unresolved_factors<'a>(value: AtomView<'a>, factors: &mut Vec<AtomView<'a>>) {
        if let AtomView::Mul(product) = value {
            for factor in product.iter() {
                Self::unresolved_factors(factor, factors);
            }
        } else if Self::has_unresolved_ports(value) {
            factors.push(value);
        }
    }

    /// Borrow unchanged branches: only bracket removal requires renormalizing
    /// their arithmetic parents.
    pub(super) fn normalize_products(value: AtomView<'_>) -> AtomOrView<'_> {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::OutputNormalization,
        );
        match value {
            AtomView::Add(sum) => {
                let terms: Vec<_> = sum.iter().map(Self::normalize_products).collect();
                if terms.iter().all(|term| matches!(term, AtomOrView::View(_))) {
                    return value.into();
                }
                // Bulk sum construction flattens nested sums, unlike the
                // existing S + S -> 2*S grouping. Retain that grouping and the
                // evaluation order of rounded coefficients and callbacks.
                if terms.iter().all(|term| {
                    !matches!(term.as_view(), AtomView::Add(_))
                        && PortRewriteState::exact_arithmetic(term.as_view())
                }) {
                    Atom::add_many(terms).into()
                } else {
                    terms
                        .iter()
                        .fold(Atom::Zero, |sum, term| sum + term.as_view())
                        .into()
                }
            }
            AtomView::Mul(product) => {
                let factors: Vec<_> = product.iter().map(Self::normalize_products).collect();
                if factors
                    .iter()
                    .all(|factor| matches!(factor, AtomOrView::View(_)))
                {
                    return value.into();
                }
                if factors
                    .iter()
                    .all(|factor| PortRewriteState::exact_arithmetic(factor.as_view()))
                {
                    Atom::mul_many(factors).into()
                } else {
                    factors
                        .iter()
                        .fold(Atom::one(), |product, factor| product * factor.as_view())
                        .into()
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                // A tensor base owns its contraction scope as a whole.
                // Removing its bracket lets scalar normalization distribute the
                // power over indexed factors and merge inner and outer dummies.
                if matches!(base, AtomView::Fun(function)
                        if function.get_symbol() == SPENSO_TAG.bracket)
                    && base.is_tensorial(StrictTensorFilter::Tagged)
                {
                    return value.into();
                }
                match Self::normalize_products(base) {
                    AtomOrView::View(_) => value.into(),
                    base => base.as_view().pow(exponent).into(),
                }
            }
            AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.bracket => {
                let mut factors = Vec::new();
                for arg in fun.iter() {
                    let factor = Self::normalize_products(arg).into_owned();
                    if let AtomView::Fun(nested) = factor.as_view()
                        && nested.get_symbol() == SPENSO_TAG.bracket
                    {
                        factors.extend(nested.iter().map(|arg| arg.to_owned()));
                    } else {
                        factors.push(factor);
                    }
                }
                let ordinary = if factors
                    .iter()
                    .all(|factor| PortRewriteState::exact_arithmetic(factor.as_view()))
                {
                    Atom::mul_many(&factors)
                } else {
                    factors
                        .iter()
                        .fold(Atom::one(), |product, factor| product * factor)
                };
                let mut before = Vec::new();
                for factor in &factors {
                    Self::unresolved_factors(factor.as_view(), &mut before);
                }
                let mut after = Vec::new();
                Self::unresolved_factors(ordinary.as_view(), &mut after);
                // Explicit indices survive sorting. Unresolved factors must
                // neither change positional order nor coalesce into powers.
                let normalized = if before != after {
                    FunctionBuilder::new(SPENSO_TAG.bracket)
                        .add_args(factors)
                        .finish()
                } else {
                    ordinary
                };
                if normalized.as_view() == value {
                    value.into()
                } else {
                    normalized.into()
                }
            }
            // Tensor arguments include scalar metadata; chain and trace arguments
            // also own placeholder scopes. Do not rewrite inside those containers.
            _ => value.into(),
        }
    }

    /// Build a presentation-only atom whose surviving tensor ports follow the
    /// public logical interface rather than `OrderedStructure`'s canonical
    /// storage order.
    pub fn presentation_atom(&self) -> Atom {
        let slots = self.structure.logical_slots();
        let mut state = PortRewriteState::new(slots.len());
        state.bulk_arithmetic = PortRewriteState::exact_arithmetic(self.expression.as_view());
        if slots.is_empty() && state.bulk_arithmetic {
            return self.expression.clone();
        }
        reorder_presentation_ports(self.expression.as_view(), &slots, &mut state, false)
    }

    /// Present an already selected fragment using this tensor's established
    /// logical interface. Additive occurrence order must survive slicing; a
    /// fragment is not a newly inferred tensor and introduces no new indices.
    pub fn presentation_fragment(&self, fragment: AtomView<'_>) -> Atom {
        let slots = self.structure.logical_slots();
        let mut state = PortRewriteState::new(slots.len());
        state.bulk_arithmetic = PortRewriteState::exact_arithmetic(fragment);
        if slots.is_empty() && state.bulk_arithmetic {
            return fragment.to_owned();
        }
        reorder_presentation_ports(fragment, &slots, &mut state, true)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PortPair {
    pub left: usize,
    pub right: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MatrixChannel {
    pub input: usize,
    pub output: usize,
}

/// Failures of port planning, admitted index history, and tensor composition.
///
/// See the [composition example](SymbolicTensor::compose).
#[derive(Debug, Error)]
pub enum TensorCompositionError {
    #[error("axes {axes:?} are not a permutation of a rank-{rank} interface")]
    InvalidAxisPermutation { rank: usize, axes: Vec<usize> },
    #[error("{0}")]
    InvalidResultInterface(String),
    #[error("{0}")]
    Inference(#[source] Box<TensorInferenceError>),
    #[error(transparent)]
    Cooking(#[from] CookingError),
    #[error(
        "outer product cannot preserve equal explicit indices on compatible port pairs {0:?}; use distinct indices"
    )]
    ExplicitOuterContraction(Vec<PortPair>),
    #[error("port position {position} is outside an interface of rank {rank}")]
    InvalidPort { position: usize, rank: usize },
    #[error("ports {left} and {right} carry incompatible representations")]
    IncompatiblePorts { left: usize, right: usize },
    #[error("explicit ports {left} and {right} do not carry the same abstract index")]
    UnequalExplicitIndices { left: usize, right: usize },
    #[error("tensor multiplication is ambiguous; compatible port pairs are {0:?}")]
    Ambiguous(Vec<PortPair>),
    #[error("the selected tensor does not have an unambiguous two-ended channel")]
    NoMatrixChannel,
    #[error("a two-ended channel requires distinct ports, got ({input}, {output})")]
    DegenerateChannel { input: usize, output: usize },
    #[error("ports ({input}, {output}) do not form an input-to-output propagation channel")]
    InvalidChannelOrientation { input: usize, output: usize },
    #[error("interface port {position} could not be located consistently in the tensor atom")]
    MissingInterfacePort { position: usize },
    #[error(
        "explicit index `{index}` occurs on more than two compatible tensor ports ({occurrences} occurrences)"
    )]
    InvalidExplicitMultiplicity {
        index: AbstractIndex,
        occurrences: usize,
    },
    #[error(transparent)]
    InvalidChainNesting(#[from] ChainNestingError),
    #[error(
        "composition would nest an existing chain or trace; only the primary channel of a root chain can be extended"
    )]
    NestedChainLike,
}

impl From<TensorInferenceError> for TensorCompositionError {
    fn from(error: TensorInferenceError) -> Self {
        match error {
            TensorInferenceError::Composition(error) => error,
            error => Self::Inference(Box::new(error)),
        }
    }
}

fn representation(slot: &PartialSlot) -> Representation<LibraryRep> {
    slot.rep()
}

fn representations_match(left: &PartialSlot, right: &PartialSlot) -> bool {
    representation(left).matches(&representation(right))
}

fn automatically_contractible(left: &PartialSlot, right: &PartialSlot) -> bool {
    if !representations_match(left, right) {
        return false;
    }

    match (left.aind, right.aind) {
        (PartialIndex::Explicit(left), PartialIndex::Explicit(right)) => left == right,
        _ => true,
    }
}

pub fn compatible_pairs(left: &PartialStructure, right: &PartialStructure) -> Vec<PortPair> {
    let left = left.logical_slots();
    let right = right.logical_slots();
    left.iter()
        .enumerate()
        .flat_map(|(left_position, left_slot)| {
            right
                .iter()
                .enumerate()
                .filter(move |(_, right_slot)| automatically_contractible(left_slot, right_slot))
                .map(move |(right_position, _)| PortPair {
                    left: left_position,
                    right: right_position,
                })
        })
        .collect()
}

fn is_structured_scalar(value: AtomView<'_>) -> bool {
    if matches!(value, AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.dot) {
        return true;
    }

    value
        .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
        .is_ok_and(|structure| structure.canonical().is_scalar())
}

fn port_matches(value: AtomView<'_>, expected: PartialSlot) -> bool {
    match expected.aind {
        PartialIndex::Explicit(index) => Slot::<LibraryRep, AbstractIndex>::try_from(value)
            .is_ok_and(|slot| slot.rep() == expected.rep() && slot.aind() == index),
        PartialIndex::Open(_) => {
            Representation::<LibraryRep>::try_from(value)
                .is_ok_and(|representation| representation == expected.rep())
                || Slot::<LibraryRep, AbstractIndex>::try_from(value).is_ok_and(|slot| {
                    slot.rep() == expected.rep()
                        && matches!(slot.aind(), AbstractIndex::Open { .. })
                })
        }
    }
}

fn chain_endpoints_match(value: AtomView<'_>, input: PartialSlot, output: PartialSlot) -> bool {
    let AtomView::Fun(function) = value else {
        return false;
    };
    if function.get_symbol() != SPENSO_TAG.chain {
        return false;
    }
    let mut arguments = function.iter();
    matches!(arguments.next(), Some(value) if port_matches(value, input))
        && matches!(arguments.next(), Some(value) if port_matches(value, output))
}

fn has_transparent_chain_channel(
    value: AtomView<'_>,
    input: PartialSlot,
    output: PartialSlot,
) -> bool {
    match value {
        AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.chain => {
            chain_endpoints_match(value, input, output)
        }
        AtomView::Fun(function) if function.get_symbol().has_tag(&SPENSO_TAG.broadcast) => {
            matches!(function.iter().collect::<Vec<_>>().as_slice(), [argument] if has_transparent_chain_channel(*argument, input, output))
        }
        AtomView::Fun(function)
            if function.get_symbol() == *shadowing::SYM
                || function.get_symbol() == *shadowing::ANTISYM
                || function.get_symbol() == *shadowing::CYCLIC =>
        {
            let tensor_arguments = function
                .iter()
                .filter(|argument| {
                    argument.is_tensorial(StrictTensorFilter::Tagged)
                        || has_transparent_chain_channel(*argument, input, output)
                })
                .collect::<Vec<_>>();
            matches!(tensor_arguments.as_slice(), [argument] if has_transparent_chain_channel(*argument, input, output))
        }
        AtomView::Mul(product) => {
            let tensor_factors = product
                .iter()
                .filter(|factor| {
                    (factor.is_tensorial(StrictTensorFilter::Tagged)
                        || matches!(
                            factor,
                            AtomView::Fun(function)
                                if function.get_symbol() == *shadowing::SYM
                                    || function.get_symbol() == *shadowing::ANTISYM
                                    || function.get_symbol() == *shadowing::CYCLIC
                        ))
                        && !is_structured_scalar(*factor)
                })
                .collect::<Vec<_>>();
            matches!(tensor_factors.as_slice(), [factor] if has_transparent_chain_channel(*factor, input, output))
        }
        AtomView::Add(sum) => {
            let terms = sum.iter().collect::<Vec<_>>();
            !terms.is_empty()
                && terms
                    .into_iter()
                    .all(|term| has_transparent_chain_channel(term, input, output))
        }
        _ => false,
    }
}

fn contains_chain_like(value: AtomView<'_>) -> bool {
    match value {
        AtomView::Add(sum) => sum.iter().any(contains_chain_like),
        AtomView::Mul(product) => product.iter().any(contains_chain_like),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            contains_chain_like(base) || contains_chain_like(exponent)
        }
        AtomView::Fun(function) => {
            function.get_symbol() == SPENSO_TAG.chain
                || function.get_symbol() == SPENSO_TAG.trace
                || function.iter().any(contains_chain_like)
        }
        _ => false,
    }
}

/// A scalar coefficient does not bind contextual ports. Preserve it while
/// extending or closing the one root chain; other wrappers retain their scope.
fn root_chain_parts(value: AtomView<'_>) -> Option<(Atom, FunView<'_>)> {
    match value {
        AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.chain => {
            Some((Atom::num(1), function))
        }
        AtomView::Mul(product) => {
            let mut chain = None;
            let mut coefficient = Atom::num(1);
            for factor in product.iter() {
                if let AtomView::Fun(function) = factor
                    && function.get_symbol() == SPENSO_TAG.chain
                    && chain.is_none()
                {
                    chain = Some(function);
                } else if !contains_chain_like(factor)
                    && (!factor.is_tensorial(StrictTensorFilter::Tagged)
                        || is_structured_scalar(factor))
                {
                    coefficient *= factor;
                } else {
                    return None;
                }
            }
            chain.map(|chain| (coefficient, chain))
        }
        _ => None,
    }
}

fn root_chain_channel_is_live(
    value: &SymbolicTensor<PartialStructure>,
    channel: MatrixChannel,
) -> bool {
    if channel
        != (MatrixChannel {
            input: 0,
            output: 1,
        })
    {
        return false;
    }
    let slots = value.structure.logical_slots();
    slots.len() >= 2
        && root_chain_parts(value.expression.as_view()).is_some_and(|(_, function)| {
            chain_endpoints_match(AtomView::Fun(function), slots[0], slots[1])
        })
}

fn validate_position(
    slots: &[PartialSlot],
    position: usize,
) -> Result<PartialSlot, TensorCompositionError> {
    slots
        .get(position)
        .copied()
        .ok_or(TensorCompositionError::InvalidPort {
            position,
            rank: slots.len(),
        })
}

fn shared_index(
    left: PartialSlot,
    right: PartialSlot,
    positions: PortPair,
) -> Result<Option<AbstractIndex>, TensorCompositionError> {
    if !representations_match(&left, &right) {
        return Err(TensorCompositionError::IncompatiblePorts {
            left: positions.left,
            right: positions.right,
        });
    }

    match (left.aind, right.aind) {
        (PartialIndex::Explicit(left), PartialIndex::Explicit(right)) if left == right => {
            Ok(Some(left))
        }
        (PartialIndex::Explicit(_), PartialIndex::Explicit(_)) => {
            Err(TensorCompositionError::UnequalExplicitIndices {
                left: positions.left,
                right: positions.right,
            })
        }
        (PartialIndex::Explicit(index), PartialIndex::Open(_))
        | (PartialIndex::Open(_), PartialIndex::Explicit(index)) => Ok(Some(index)),
        (PartialIndex::Open(_), PartialIndex::Open(_)) => Ok(None),
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Reserve all written and established interface indices once before a
    /// sequence of allocations. The parser owns serialized collision handling.
    pub fn reserved_dummies<'a>(
        values: impl IntoIterator<Item = &'a Self>,
    ) -> ParseState<AbstractIndex> {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::DummyReservation,
        );
        let state = ParseState::default();
        for value in values {
            value.proofs.reserve_indices(
                &value.expression,
                value
                    .structure
                    .logical_slots()
                    .into_iter()
                    .filter_map(|slot| match slot.aind {
                        PartialIndex::Explicit(index) => Some(index),
                        PartialIndex::Open(_) => None,
                    }),
                &state,
            );
        }
        state
    }
}

impl super::TensorProofs {
    /// Reserve source names through the same borrowed observer for either interface.
    pub(super) fn reserve_indices(
        &self,
        expression: &Atom,
        explicit_indices: impl IntoIterator<Item = AbstractIndex>,
        state: &ParseState<AbstractIndex>,
    ) {
        let explicit_indices = explicit_indices.into_iter().collect::<Vec<_>>();
        let public_indices = explicit_indices.iter().copied().collect::<HashSet<_>>();
        if let Some(occurrences) = self.occurrences.get() {
            for (slot, _) in occurrences.entries() {
                if !public_indices.contains(&slot.aind()) {
                    state.reserve_index(slot.aind());
                }
            }
        }
        if let Some(observed) = self.observations.get() {
            for &index in &observed.reserved_indices {
                state.reserve_index(index);
            }
        } else {
            #[cfg(test)]
            DUMMY_RESERVATION_WALKS.with(|count| count.set(count.get() + 1));
            state.reserve_indices(expression.as_view());
        }
        for index in explicit_indices {
            state.reserve_index(index);
        }
    }
}

pub(super) fn is_composite_head(symbol: Symbol, tags: &SpensoTags) -> bool {
    symbol == tags.bracket
        || symbol == tags.chain
        || symbol == tags.trace
        || symbol.has_tag(&tags.broadcast)
        // Registered projectors carry these attributes. Check the bits before
        // accessing their lazy bundles, which each synchronize initialization.
        || (symbol.is_symmetric() && symbol == *shadowing::SYM)
        || (symbol.is_antisymmetric() && symbol == *shadowing::ANTISYM)
        || (symbol.is_cyclesymmetric() && symbol == *shadowing::CYCLIC)
}

pub(super) fn is_tensor_leaf_head(symbol: Symbol) -> bool {
    let tags = &*SPENSO_TAG;
    symbol.has_tag(&tags.tensor) && !is_composite_head(symbol, tags)
}

fn direct_structural_port(value: AtomView<'_>) -> bool {
    Slot::<LibraryRep, AbstractIndex>::try_from(value).is_ok()
        || Representation::<LibraryRep>::try_from(value).is_ok()
}

/// Explicit identities visible in one multiplicative scope, in total key order.
///
/// The aligned mask owns collision history: an unset bit means one occurrence,
/// a set bit means two uses, which may be two free ports. Cancellation and compact
/// notation do not erase it. Opposite orientations share a base key, while an
/// entry retains its original slot orientation (exact for a singleton). The free
/// structure separately owns port roles, orientation and logical layout.
/// Additive and closed power scopes hide their
/// internal identities before exporting their boundary observations.
#[derive(Clone, Debug)]
pub(super) struct ExplicitIndexOccurrences {
    indices: Vec<Slot<LibraryRep, AbstractIndex>>,
    dummies: SubSet<SlotIndex>,
    /// Per-branch base/dual maxima, packed into two bits per orientation.
    /// Alternatives can have both maxima present without adding their totals.
    orientations: Vec<u8>,
    /// Admission defers its first third-occurrence error until scope validation
    /// finishes. Invalid local branches must not disappear during scoping.
    overflow: Option<Slot<LibraryRep, AbstractIndex>>,
    /// A queried pair may live in a closed power scope rather than the outer scope.
    pub(super) selected_pair: bool,
}

#[derive(Clone, Copy, Debug)]
enum OccurrenceMerge {
    Product,
    Alternatives,
}

impl Default for ExplicitIndexOccurrences {
    fn default() -> Self {
        Self {
            indices: Vec::new(),
            dummies: iter::empty().collect(),
            orientations: Vec::new(),
            overflow: None,
            selected_pair: false,
        }
    }
}

impl ExplicitIndexOccurrences {
    pub(super) fn from_atom(value: AtomView<'_>, slots: &mut SlotMatcher) -> Self {
        InterfaceInference::index_occurrences(value, slots, &mut AHashMap::new())
    }

    /// Combine admitted construction metadata without inspecting result syntax.
    pub(super) fn product(
        sources: &[&SymbolicTensor<PartialStructure>],
    ) -> Result<Self, TensorCompositionError> {
        let mut result = Self::default();
        for source in sources {
            let positional = !source.structure.open_positions().is_empty();
            let raw_observation;
            let observed = if positional {
                if source.expression.is_zero() && source.proofs.occurrences.get().is_none() {
                    let boundary = Self::from_interface(&source.structure);
                    boundary.validate()?;
                    let _ = source.proofs.occurrences.set(Arc::new(boundary));
                }
                if let Some(observed) = source.proofs.occurrences.get() {
                    observed.as_ref()
                } else {
                    // A raw positional layout may differ from physical argument
                    // order. Observe scopes locally; this is not authoritative
                    // metadata and must never populate the admitted proof cell.
                    raw_observation = InterfaceInference::scoped_occurrences(
                        source.expression.as_view(),
                        &mut SlotMatcher::default(),
                    )
                    .map_err(TensorCompositionError::from)?;
                    raw_observation.validate()?;
                    &raw_observation
                }
            } else {
                source
                    .explicit_index_occurrences()
                    .map_err(TensorCompositionError::from)?
            };
            if positional {
                let mut explicit = Self {
                    selected_pair: observed.selected_pair,
                    ..Self::default()
                };
                for (slot, count, base, dual) in observed.oriented_entries() {
                    if !CookSettings::index_has_open_port(slot.aind())
                        .map_err(TensorCompositionError::from)?
                    {
                        explicit.add_oriented(slot, count, base, dual);
                    }
                }
                result.append(&explicit);
            } else {
                result.append(observed);
            }
        }
        result.validate()?;
        Ok(result)
    }

    /// Seed explicitly declared public ports, including the nominal ports of zero.
    pub(super) fn from_interface(interface: &PartialStructure) -> Self {
        let mut result = Self::default();
        for slot in interface.logical_slots() {
            if let PartialIndex::Explicit(index) = slot.aind {
                result.add(slot.rep().slot(index), 1);
            }
        }
        result
    }

    /// Export only identities which can meet the surrounding multiplicative scope.
    pub(super) fn scoped(&self, interface: &PartialStructure) -> Self {
        let boundary = Self::from_interface(interface);
        let mut result = self.clone();
        result.retain(|slot| boundary.count(*slot) != 0);
        result.branch_max(&boundary);
        result
    }

    /// Transport every retained identity with the same mapping as the free ports.
    pub(super) fn mapped<E>(
        &self,
        mut map: impl FnMut(
            Slot<LibraryRep, AbstractIndex>,
        ) -> Result<Slot<LibraryRep, AbstractIndex>, E>,
    ) -> Result<Self, E> {
        let mut result = Self {
            selected_pair: self.selected_pair,
            ..Self::default()
        };
        for (slot, total, base, dual) in self.oriented_entries() {
            let representation = slot.rep().base();
            let mapped_base = if base != 0 {
                Some(map(representation.slot(slot.aind()))?)
            } else {
                None
            };
            let mapped_dual = if dual != 0 {
                Some(map(representation.dual().slot(slot.aind()))?)
            } else {
                None
            };
            match (mapped_base, mapped_dual) {
                (Some(left), Some(right)) if Self::compare(&left, &right).is_eq() => {
                    let mut bounds = [0usize; 2];
                    for (slot, count) in [(left, base), (right, dual)] {
                        let orientation = usize::from(!slot.rep().rep.is_base());
                        bounds[orientation] = bounds[orientation].saturating_add(count).min(total);
                    }
                    let representative = if !slot.rep().rep.is_base() {
                        right
                    } else {
                        left
                    };
                    result.add_oriented(representative, total, bounds[0], bounds[1]);
                }
                (left, right) => {
                    if let Some(slot) = left {
                        result.add(slot, base);
                    }
                    if let Some(slot) = right {
                        result.add(slot, dual);
                    }
                }
            }
        }
        Ok(result)
    }

    /// Dimension's structural ordering deliberately ties symbolic dimensions.
    /// Collision keys instead need equality to be equivalent to comparator ties.
    fn compare(
        left: &Slot<LibraryRep, AbstractIndex>,
        right: &Slot<LibraryRep, AbstractIndex>,
    ) -> CmpOrdering {
        left.rep()
            .base()
            .rep
            .cmp(&right.rep().base().rep)
            .then_with(|| match (left.rep().dim, right.rep().dim) {
                (Dimension::Concrete(left), Dimension::Concrete(right)) => left.cmp(&right),
                (Dimension::Concrete(_), Dimension::Symbolic(_)) => CmpOrdering::Less,
                (Dimension::Symbolic(_), Dimension::Concrete(_)) => CmpOrdering::Greater,
                (Dimension::Symbolic(left), Dimension::Symbolic(right)) => left.cmp(&right),
            })
            .then_with(|| left.aind().cmp(&right.aind()))
    }

    pub(super) fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }

    pub(super) fn entries(
        &self,
    ) -> impl Iterator<Item = (Slot<LibraryRep, AbstractIndex>, usize)> + '_ {
        self.indices
            .iter()
            .copied()
            .zip(self.dummies.iter())
            .map(|(slot, dummy)| {
                let count = if self
                    .overflow
                    .is_some_and(|overflow| Self::compare(&overflow, &slot).is_eq())
                {
                    3
                } else {
                    1 + usize::from(dummy)
                };
                (slot, count)
            })
    }

    /// Total multiplicity and independent per-variance branch bounds.
    pub(super) fn oriented_entries(
        &self,
    ) -> impl Iterator<Item = (Slot<LibraryRep, AbstractIndex>, usize, usize, usize)> + '_ {
        self.entries()
            .zip(&self.orientations)
            .map(|((slot, total), &orientations)| {
                (
                    slot,
                    total,
                    usize::from(orientations & 3),
                    usize::from((orientations >> 2) & 3),
                )
            })
    }

    /// A replacement may remove internal pairs, but cannot introduce a new
    /// encoded index or increase its maximum count in an additive branch.
    pub(super) fn is_bounded_by(&self, source: &Self) -> bool {
        self.oriented_entries().all(|(slot, count, base, dual)| {
            let Ok(position) = source
                .indices
                .binary_search_by(|candidate| Self::compare(candidate, &slot))
            else {
                return false;
            };
            count <= source.count(slot)
                && base <= usize::from(source.orientations[position] & 3)
                && dual <= usize::from((source.orientations[position] >> 2) & 3)
        })
    }

    pub(super) fn count(&self, target: Slot<LibraryRep, AbstractIndex>) -> usize {
        let target = target.rep().base().slot(target.aind());
        match self
            .indices
            .binary_search_by(|slot| Self::compare(slot, &target))
        {
            Ok(_)
                if self
                    .overflow
                    .is_some_and(|overflow| Self::compare(&overflow, &target).is_eq()) =>
            {
                3
            }
            Ok(position) => 1 + usize::from(self.dummies[SlotIndex(position)]),
            Err(_) => 0,
        }
    }

    pub(super) fn add(&mut self, slot: Slot<LibraryRep, AbstractIndex>, count: usize) {
        let (base, dual) = if !slot.rep().rep.is_base() {
            (0, count)
        } else {
            (count, 0)
        };
        self.add_oriented(slot, count, base, dual);
    }

    pub(super) fn add_oriented(
        &mut self,
        slot: Slot<LibraryRep, AbstractIndex>,
        count: usize,
        base: usize,
        dual: usize,
    ) {
        if count == 0 {
            return;
        }
        match self
            .indices
            .binary_search_by(|candidate| Self::compare(candidate, &slot))
        {
            Ok(position) => {
                let previous = 1 + usize::from(self.dummies[SlotIndex(position)]);
                if previous.saturating_add(count) > 2 {
                    self.overflow.get_or_insert(slot);
                }
                self.dummies.add(SlotIndex(position));
                let orientations = self.orientations[position];
                self.orientations[position] =
                    (usize::from(orientations & 3).saturating_add(base).min(3) as u8)
                        | ((usize::from((orientations >> 2) & 3)
                            .saturating_add(dual)
                            .min(3) as u8)
                            << 2);
            }
            Err(position) => {
                let append = position == self.indices.len();
                self.indices.insert(position, slot);
                self.orientations
                    .insert(position, (base.min(3) as u8) | ((dual.min(3) as u8) << 2));
                if append {
                    self.dummies.push(count >= 2);
                } else {
                    self.dummies = self
                        .dummies
                        .iter()
                        .take(position)
                        .chain(iter::once(count >= 2))
                        .chain(self.dummies.iter().skip(position))
                        .collect();
                }
                if count > 2 {
                    self.overflow.get_or_insert(slot);
                }
            }
        }
    }

    pub(super) fn append(&mut self, other: &Self) {
        self.merge(other, OccurrenceMerge::Product);
    }

    /// Remove one selected public occurrence before simultaneous relabelling.
    fn remove(
        &mut self,
        slot: Slot<LibraryRep, AbstractIndex>,
    ) -> Result<(), TensorCompositionError> {
        self.validate()?;
        let target = slot.rep().base().slot(slot.aind());
        let position = self
            .indices
            .binary_search_by(|candidate| Self::compare(candidate, &target))
            .map_err(|_| {
                TensorCompositionError::InvalidResultInterface(
                    "selected explicit port is absent from admitted index metadata".into(),
                )
            })?;
        if self.dummies[SlotIndex(position)] {
            let orientation = u32::from(!slot.rep().rep.is_base()) * 2;
            let mut bounds = self.orientations[position];
            let count = (bounds >> orientation) & 3;
            if count == 0 {
                return Err(TensorCompositionError::InvalidResultInterface(
                    "selected explicit port orientation is absent from admitted index metadata"
                        .into(),
                ));
            }
            bounds = (bounds & !(3 << orientation)) | ((count - 1) << orientation);
            self.orientations[position] = (bounds & 3).min(1) | ((((bounds >> 2) & 3).min(1)) << 2);
            let stored_orientation = u32::from(!self.indices[position].rep().rep.is_base()) * 2;
            if (self.orientations[position] >> stored_orientation) & 3 == 0 {
                let old = self.indices[position];
                let representation = old.rep().base();
                self.indices[position] = if self.orientations[position] & 3 != 0 {
                    representation.slot(old.aind())
                } else {
                    representation.dual().slot(old.aind())
                };
            }
            self.dummies.sub(SlotIndex(position));
        } else {
            self.indices.remove(position);
            self.orientations.remove(position);
            self.dummies = self
                .dummies
                .iter()
                .enumerate()
                .filter_map(|(index, dummy)| (index != position).then_some(dummy))
                .collect();
        }
        Ok(())
    }

    pub(super) fn branch_max(&mut self, other: &Self) {
        self.merge(other, OccurrenceMerge::Alternatives);
    }

    fn merge(&mut self, other: &Self, mode: OccurrenceMerge) {
        self.selected_pair |= other.selected_pair;
        let mut indices = Vec::with_capacity(self.indices.len() + other.indices.len());
        let mut dummies: SubSet<SlotIndex> = iter::empty().collect();
        let mut orientations = Vec::with_capacity(self.indices.len() + other.indices.len());
        let mut overflow = self.overflow.or(other.overflow);
        let mut left = self.oriented_entries().peekable();
        let mut right = other.oriented_entries().peekable();
        loop {
            let (slot, count, base, dual) = match (left.peek(), right.peek()) {
                (Some(&(left_slot, _, _, _)), Some(&(right_slot, _, _, _))) => {
                    match Self::compare(&left_slot, &right_slot) {
                        CmpOrdering::Less => left
                            .next()
                            .expect("peek established a remaining left index"),
                        CmpOrdering::Greater => right
                            .next()
                            .expect("peek established a remaining right index"),
                        CmpOrdering::Equal => {
                            let (slot, left_count, left_base, left_dual) = left
                                .next()
                                .expect("peek established a remaining left index");
                            let (_, right_count, right_base, right_dual) = right
                                .next()
                                .expect("peek established a remaining right index");
                            (
                                slot,
                                if matches!(mode, OccurrenceMerge::Alternatives) {
                                    left_count.max(right_count)
                                } else {
                                    left_count.saturating_add(right_count)
                                },
                                if matches!(mode, OccurrenceMerge::Alternatives) {
                                    left_base.max(right_base)
                                } else {
                                    left_base.saturating_add(right_base)
                                },
                                if matches!(mode, OccurrenceMerge::Alternatives) {
                                    left_dual.max(right_dual)
                                } else {
                                    left_dual.saturating_add(right_dual)
                                },
                            )
                        }
                    }
                }
                (Some(_), None) => left
                    .next()
                    .expect("peek established a remaining left index"),
                (None, Some(_)) => right
                    .next()
                    .expect("peek established a remaining right index"),
                (None, None) => break,
            };
            indices.push(slot);
            dummies.push(count >= 2);
            orientations.push((base.min(3) as u8) | ((dual.min(3) as u8) << 2));
            if count > 2 {
                overflow.get_or_insert(slot);
            }
        }
        drop((left, right));
        self.indices = indices;
        self.dummies = dummies;
        self.orientations = orientations;
        self.overflow = overflow;
    }

    pub(super) fn retain(
        &mut self,
        mut keep: impl FnMut(&Slot<LibraryRep, AbstractIndex>) -> bool,
    ) {
        let mut position = 0;
        let mut dummies: SubSet<SlotIndex> = iter::empty().collect();
        let mut orientations = Vec::new();
        self.indices.retain(|slot| {
            let retain = keep(slot);
            if retain {
                dummies.push(self.dummies[SlotIndex(position)]);
                orientations.push(self.orientations[position]);
            }
            position += 1;
            retain
        });
        self.dummies = dummies;
        self.orientations = orientations;
    }

    pub(super) fn repeat(&mut self, repetitions: usize) {
        if repetitions == 0 {
            self.indices.clear();
            self.dummies = iter::empty().collect();
            self.orientations.clear();
        } else if repetitions > 1 {
            self.overflow = self.overflow.or_else(|| {
                self.entries()
                    .find(|(_, count)| count.saturating_mul(repetitions) > 2)
                    .map(|(slot, _)| slot)
            });
            self.dummies = SubSet::full(self.indices.len());
            for bounds in &mut self.orientations {
                *bounds = (usize::from(*bounds & 3).saturating_mul(repetitions).min(3) as u8)
                    | ((usize::from((*bounds >> 2) & 3)
                        .saturating_mul(repetitions)
                        .min(3) as u8)
                        << 2);
            }
        }
    }

    pub(super) fn validate(&self) -> Result<(), TensorCompositionError> {
        if let Some(slot) = self.overflow {
            return Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: slot.aind(),
                occurrences: 3,
            });
        }
        Ok(())
    }
}

/// Reject Einstein indices that occur more than twice in any additive branch.
///
/// The walk follows tensor-bearing arguments only, so representation-valued
/// scalar metadata is not mistaken for a structural port. Additive branches
/// are checked independently because each summand owns the same public
/// interface rather than contributing another occurrence of it.
pub fn validate_explicit_index_occurrences(atom: &Atom) -> Result<(), TensorCompositionError> {
    let occurrences =
        ExplicitIndexOccurrences::from_atom(atom.as_view(), &mut SlotMatcher::default());
    occurrences.validate()
}

#[derive(Clone)]
struct PortRewriteState {
    claimed: Vec<bool>,
    applied: Vec<bool>,
    bulk_arithmetic: bool,
    /// Admitted explicit metadata can expose copies compressed into one power.
    explicit_copies: bool,
    indices: Option<ParseState<AbstractIndex>>,
}

impl PortRewriteState {
    fn new(rank: usize) -> Self {
        Self {
            claimed: vec![false; rank],
            applied: vec![false; rank],
            bulk_arithmetic: false,
            explicit_copies: false,
            indices: None,
        }
    }

    /// Regroup exact arithmetic only when rebuilding cannot invoke user code.
    /// Metric normalization is intrinsic; callbacks in its operands are checked too.
    fn exact_arithmetic(value: AtomView<'_>) -> bool {
        let mut exact = true;
        let mut slots = SlotMatcher::default();
        value.visitor(&mut |node| {
            exact &= !node.needs_normalization()
                && match node {
                    AtomView::Num(number) => matches!(
                        number.get_coeff_view(),
                        CoefficientView::Natural(..) | CoefficientView::Large(..)
                    ),
                    AtomView::Fun(function) => {
                        let symbol = function.get_symbol();
                        symbol.get_evaluation_info().is_none()
                            && (symbol == ETS.metric
                                || symbol.get_normalization_function().is_none())
                    }
                    _ => true,
                };
            if exact && let SlotMatch::Explicit(slot) = slots.classify(node) {
                // Relabelling replaces the complete admitted port. Its index
                // payload is an opaque identity and is never rebuilt here.
                exact = Self::exact_arithmetic(slot.dimension());
                return false;
            }
            exact
        });
        exact
    }
}

pub(crate) fn matching_interface_position(
    value: AtomView<'_>,
    slots: &[PartialSlot],
    claimed: &[bool],
) -> Option<usize> {
    if let Ok(slot) = Slot::<LibraryRep, AbstractIndex>::try_from(value) {
        return slots.iter().enumerate().position(|(position, expected)| {
            if claimed[position] || expected.rep() != slot.rep() {
                return false;
            }
            match expected.aind {
                PartialIndex::Explicit(index) => index == slot.aind(),
                // Atom-level open markers are local to a tensor occurrence,
                // while partial-interface IDs are renumbered whenever
                // interfaces are combined or ports are removed.
                PartialIndex::Open(_) => matches!(slot.aind(), AbstractIndex::Open { .. }),
            }
        });
    }

    let representation = Representation::<LibraryRep>::try_from(value).ok()?;
    slots.iter().enumerate().position(|(position, expected)| {
        !claimed[position]
            && expected.rep() == representation
            && matches!(expected.aind, PartialIndex::Open(_))
    })
}

fn reorder_presentation_ports(
    value: AtomView<'_>,
    slots: &[PartialSlot],
    state: &mut PortRewriteState,
    preserve_additive_order: bool,
) -> Atom {
    match value {
        AtomView::Add(add) => {
            let initial = state.clone();
            let mut common = None::<PortRewriteState>;
            let terms: Vec<_> = add
                .iter()
                .map(|term| {
                    let mut local = initial.clone();
                    let term = reorder_presentation_ports(term, slots, &mut local, true);
                    if let Some(common) = &mut common {
                        for (common, local) in common.claimed.iter_mut().zip(&local.claimed) {
                            *common &= local;
                        }
                    } else {
                        common = Some(local);
                    }
                    term
                })
                .collect();
            *state = common.unwrap_or(initial);
            if state.bulk_arithmetic && terms.iter().zip(add.iter()).all(|(a, b)| a.as_view() == b)
            {
                value.to_owned()
            } else if state.bulk_arithmetic
                && terms
                    .iter()
                    .all(|term| !matches!(term.as_view(), AtomView::Add(_)))
            {
                Atom::add_many(terms)
            } else {
                terms.into_iter().fold(Atom::Zero, |sum, term| sum + term)
            }
        }
        AtomView::Mul(mul) => {
            let factors: Vec<_> = mul
                .iter()
                .map(|factor| {
                    let tensorial = factor.is_tensorial(StrictTensorFilter::Tagged)
                        || matches!(
                            factor,
                            AtomView::Fun(fun)
                                if fun.get_symbol() == *shadowing::SYM
                                    || fun.get_symbol() == *shadowing::ANTISYM
                                    || fun.get_symbol() == *shadowing::CYCLIC
                        );
                    if tensorial {
                        reorder_presentation_ports(factor, slots, state, preserve_additive_order)
                    } else {
                        factor.to_owned()
                    }
                })
                .collect();
            if state.bulk_arithmetic
                && factors
                    .iter()
                    .zip(mul.iter())
                    .all(|(a, b)| a.as_view() == b)
            {
                value.to_owned()
            } else if state.bulk_arithmetic {
                Atom::mul_many(factors)
            } else {
                factors
                    .into_iter()
                    .fold(Atom::one(), |product, factor| product * factor)
            }
        }
        AtomView::Pow(pow) => {
            let (base, exponent) = pow.get_base_exp();
            let rewritten = reorder_presentation_ports(base, slots, state, preserve_additive_order);
            if state.bulk_arithmetic && rewritten.as_view() == base {
                value.to_owned()
            } else {
                rewritten.pow(exponent.to_owned())
            }
        }
        // A compact dot has no public ports of its own.
        AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.dot => value.to_owned(),
        AtomView::Fun(fun) if is_tensor_leaf_head(fun.get_symbol()) => {
            let mut arguments = fun
                .iter()
                .map(|argument| argument.to_owned())
                .collect::<Vec<_>>();
            let matched = arguments
                .iter()
                .enumerate()
                .filter_map(|(argument_position, argument)| {
                    let interface_position =
                        matching_interface_position(argument.as_view(), slots, &state.claimed)?;
                    state.claimed[interface_position] = true;
                    Some((argument_position, interface_position, argument.clone()))
                })
                .collect::<Vec<_>>();
            // Within a sum, each term's explicit argument order records its map to the
            // shared interface; normalizing it would turn A(i,j) + A(j,i) into A(i,j) + A(i,j).
            if !preserve_additive_order {
                let mut ordered = matched.clone();
                ordered.sort_by_key(|(_, interface_position, _)| *interface_position);
                for ((argument_position, _, _), (_, _, argument)) in
                    matched.into_iter().zip(ordered)
                {
                    arguments[argument_position] = argument;
                }
            }
            if state.bulk_arithmetic
                && arguments
                    .iter()
                    .zip(fun.iter())
                    .all(|(a, b)| a.as_view() == b)
            {
                value.to_owned()
            } else {
                FunctionBuilder::new(fun.get_symbol())
                    .add_args(arguments)
                    .finish()
            }
        }
        AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.trace => {
            let mut arguments = fun.iter();
            let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
            if let Some(representation) = arguments.next() {
                rebuilt = rebuilt.add_arg(representation);
            }
            for argument in arguments {
                rebuilt = rebuilt.add_arg(reorder_presentation_ports(
                    argument,
                    slots,
                    state,
                    preserve_additive_order,
                ));
            }
            rebuilt.finish()
        }
        AtomView::Fun(fun) => {
            let symbol = fun.get_symbol();
            let mut rebuilt = FunctionBuilder::new(symbol);
            for argument in fun.iter() {
                if direct_structural_port(argument) {
                    if let Some(position) =
                        matching_interface_position(argument, slots, &state.claimed)
                    {
                        state.claimed[position] = true;
                    }
                    rebuilt = rebuilt.add_arg(argument);
                } else if is_composite_head(symbol, &SPENSO_TAG)
                    && (argument.is_tensorial(StrictTensorFilter::Tagged)
                        || matches!(
                            argument,
                            AtomView::Fun(nested)
                                if nested.get_symbol() == *shadowing::SYM
                                    || nested.get_symbol() == *shadowing::ANTISYM
                                    || nested.get_symbol() == *shadowing::CYCLIC
                        ))
                {
                    rebuilt = rebuilt.add_arg(reorder_presentation_ports(
                        argument,
                        slots,
                        state,
                        preserve_additive_order,
                    ));
                } else {
                    rebuilt = rebuilt.add_arg(argument);
                }
            }
            rebuilt.finish()
        }
        _ => value.to_owned(),
    }
}

fn rewrite_ports<'a>(
    value: AtomView<'a>,
    slots: &[PartialSlot],
    replacements: &HashMap<usize, Atom>,
    state: &mut PortRewriteState,
) -> AtomOrView<'a> {
    match value {
        AtomView::Add(add) => {
            let original = state.clone();
            let mut initial = original.clone();
            let scoped = if state.explicit_copies
                && !value.needs_normalization()
                && !replacements.is_empty()
            {
                InterfaceInference::replacement_interface(value)
                    .ok()
                    .map(|interface| {
                        let mut allowed = vec![false; slots.len()];
                        let mut claimed = original.claimed.clone();
                        for slot in interface.logical_slots() {
                            if let Some(position) =
                                slots.iter().enumerate().position(|(position, candidate)| {
                                    !claimed[position] && *candidate == slot
                                })
                            {
                                claimed[position] = true;
                                allowed[position] = true;
                            }
                        }
                        allowed
                    })
            } else {
                None
            };
            if let Some(allowed) = &scoped {
                if !allowed
                    .iter()
                    .enumerate()
                    .any(|(position, &allowed)| allowed && replacements.contains_key(&position))
                {
                    return value.into();
                }
                for (position, claimed) in initial.claimed.iter_mut().enumerate() {
                    *claimed |= !allowed[position];
                }
            }
            let mut common = None::<PortRewriteState>;
            let mut sum = value.needs_normalization().then_some(Atom::Zero);
            let mut terms = state.bulk_arithmetic.then(Vec::new);
            for (position, term) in add.iter().enumerate() {
                let mut local = initial.clone();
                let term = rewrite_ports(term, slots, replacements, &mut local);
                if let Some(common) = &mut common {
                    for (common, local) in common.claimed.iter_mut().zip(&local.claimed) {
                        *common &= local;
                    }
                    for (common, local) in common.applied.iter_mut().zip(&local.applied) {
                        *common &= local;
                    }
                } else {
                    common = Some(local);
                }
                if let Some(terms) = &mut terms {
                    terms.push(term);
                } else if sum.is_some() || !matches!(term, AtomOrView::View(_)) {
                    // Earlier borrowed terms are byte-identical and had no
                    // callbacks. Preserve the old fold order once anything changes.
                    let prefix = sum.unwrap_or_else(|| {
                        add.iter()
                            .take(position)
                            .fold(Atom::Zero, |sum, term| sum + term)
                    });
                    sum = Some(prefix + term.as_view());
                }
            }
            *state = common.unwrap_or(initial);
            if let Some(allowed) = scoped {
                for (position, allowed) in allowed.into_iter().enumerate() {
                    if !allowed {
                        state.claimed[position] = original.claimed[position];
                        state.applied[position] = original.applied[position];
                    }
                }
            }
            if let Some(terms) = terms {
                return if terms.iter().all(|term| matches!(term, AtomOrView::View(_))) {
                    value.into()
                } else if !terms.iter().all(|term| {
                    !matches!(term.as_view(), AtomView::Add(_))
                        && PortRewriteState::exact_arithmetic(term.as_view())
                }) {
                    terms
                        .iter()
                        .fold(Atom::Zero, |sum, term| sum + term.as_view())
                        .into()
                } else {
                    Atom::add_many(terms).into()
                };
            }
            sum.map_or_else(|| value.into(), Into::into)
        }
        AtomView::Mul(mul) => {
            let mut product = value.needs_normalization().then(|| Atom::num(1));
            let mut factors = state.bulk_arithmetic.then(Vec::new);
            for (position, factor) in mul.iter().enumerate() {
                let tensorial = factor.is_tensorial(StrictTensorFilter::Tagged)
                    || matches!(
                        factor,
                        AtomView::Fun(fun)
                            if fun.get_symbol() == *SYM
                                || fun.get_symbol() == *ANTISYM
                                || fun.get_symbol() == *CYCLIC
                    );
                let factor = if tensorial {
                    rewrite_ports(factor, slots, replacements, state)
                } else {
                    factor.into()
                };
                if let Some(factors) = &mut factors {
                    factors.push(factor);
                } else if product.is_some() || !matches!(factor, AtomOrView::View(_)) {
                    let prefix = product.unwrap_or_else(|| {
                        mul.iter()
                            .take(position)
                            .fold(Atom::num(1), |product, factor| product * factor)
                    });
                    product = Some(prefix * factor.as_view());
                }
            }
            if let Some(factors) = factors {
                return if factors
                    .iter()
                    .all(|factor| matches!(factor, AtomOrView::View(_)))
                {
                    value.into()
                } else if !factors
                    .iter()
                    .all(|factor| PortRewriteState::exact_arithmetic(factor.as_view()))
                {
                    factors
                        .iter()
                        .fold(Atom::one(), |product, factor| product * factor.as_view())
                        .into()
                } else {
                    Atom::mul_many(factors).into()
                };
            }
            product.map_or_else(|| value.into(), Into::into)
        }
        AtomView::Pow(pow) => {
            let (base, exponent) = pow.get_base_exp();
            let repetitions = Rational::try_from(exponent)
                .ok()
                .filter(|exponent| exponent.denominator() == 1)
                .and_then(|exponent| usize::try_from(exponent.numerator().clone()).ok())
                .filter(|&count| count > 1 && count <= slots.len());
            let observed = if state.explicit_copies
                && !value.needs_normalization()
                && !replacements.is_empty()
            {
                InterfaceInference::replacement_interface(base).ok()
            } else {
                None
            };
            if observed
                .as_ref()
                .is_some_and(|interface| interface.canonical().is_scalar())
            {
                // A closed base owns its internal identities; matching labels
                // outside this power do not make them public rewrite targets.
                return value.into();
            }
            if state.explicit_copies
                && !value.needs_normalization()
                && !replacements.is_empty()
                && let Some(repetitions) = repetitions
                && let Some(interface) = observed
            {
                let boundary = interface.logical_slots();
                if !boundary.is_empty()
                    && boundary.iter().all(|slot| match slot.aind {
                        PartialIndex::Explicit(index) => {
                            CookSettings::index_has_open_port(index).is_ok_and(|open| !open)
                        }
                        PartialIndex::Open(_) => false,
                    })
                    && boundary
                        .len()
                        .checked_mul(repetitions)
                        .is_some_and(|count| count <= slots.len())
                {
                    // A normalized power may compress several public copies.
                    // Plan complete base interfaces in copy order, reserving
                    // unmodified axes too, before invoking any port callback.
                    let mut claimed = state.claimed.clone();
                    let copies = (0..repetitions)
                        .map(|_| {
                            boundary
                                .iter()
                                .map(|slot| {
                                    let position = slots.iter().enumerate().position(
                                        |(position, candidate)| {
                                            !claimed[position] && candidate == slot
                                        },
                                    )?;
                                    claimed[position] = true;
                                    Some(position)
                                })
                                .collect::<Option<Vec<_>>>()
                        })
                        .collect::<Option<Vec<_>>>();
                    if let Some(copies) = copies
                        && copies
                            .iter()
                            .flatten()
                            .any(|position| replacements.contains_key(position))
                    {
                        let mut factors = Vec::with_capacity(repetitions);
                        for positions in copies {
                            let mut local = state.clone();
                            for (position, claimed) in local.claimed.iter_mut().enumerate() {
                                *claimed = !positions.contains(&position);
                            }
                            local.indices = state.indices.take();
                            factors.push(rewrite_ports(base, slots, replacements, &mut local));
                            for position in positions {
                                state.claimed[position] = local.claimed[position];
                                state.applied[position] |= local.applied[position];
                            }
                            state.indices = local.indices;
                        }
                        return if factors
                            .iter()
                            .all(|factor| matches!(factor, AtomOrView::View(_)))
                        {
                            value.into()
                        } else if state.bulk_arithmetic
                            && factors
                                .iter()
                                .all(|factor| PortRewriteState::exact_arithmetic(factor.as_view()))
                        {
                            Atom::mul_many(factors).into()
                        } else {
                            factors
                                .iter()
                                .fold(Atom::one(), |product, factor| product * factor.as_view())
                                .into()
                        };
                    }
                }
            }
            let base = rewrite_ports(base, slots, replacements, state);
            if !value.needs_normalization() && matches!(base, AtomOrView::View(_)) {
                value.into()
            } else {
                base.as_view().pow(exponent).into()
            }
        }
        AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.trace => {
            let mut args = fun.iter();
            let mut rewritten = Vec::with_capacity(fun.get_nargs());
            if let Some(rep) = args.next() {
                rewritten.push(rep.into());
            }
            for arg in args {
                rewritten.push(rewrite_ports(arg, slots, replacements, state));
            }
            if !value.needs_normalization()
                && fun.get_symbol().get_normalization_function().is_none()
                && fun.get_symbol().get_evaluation_info().is_none()
                && rewritten
                    .iter()
                    .all(|arg| matches!(arg, AtomOrView::View(_)))
            {
                value.into()
            } else {
                FunctionBuilder::new(fun.get_symbol())
                    .add_args(rewritten)
                    .finish()
                    .into()
            }
        }
        // A dot's rank-one channel is internal and therefore does not consume
        // positions in the dot expression's external interface.
        AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.dot => value.into(),
        AtomView::Fun(fun) => {
            if fun.get_symbol() == SPENSO_TAG.chain
                && fun.iter().take(2).any(|endpoint| {
                    matching_interface_position(endpoint, slots, &state.claimed)
                        .and_then(|position| replacements.get(&position))
                        .is_some_and(|replacement| {
                            matches!(replacement.as_view(), AtomView::Fun(vector)
                            if vector.get_symbol().has_tag(&SPENSO_TAG.rank1))
                        })
                })
            {
                // Apply substitutions through the indexed matrix factors.
                // Structural contraction can recollect their supplied endpoints.
                let indices = state.indices.get_or_insert_with(ParseState::default);
                indices.reserve_indices(value);
                if let Ok(indexed) = indices.materialize_indexed_chain(fun) {
                    return rewrite_ports(indexed.as_view(), slots, replacements, state)
                        .into_owned()
                        .into();
                }
            }
            let tensor_leaf = is_tensor_leaf_head(fun.get_symbol());
            let mut rewritten = Vec::with_capacity(fun.get_nargs());
            for (position, arg) in fun.iter().enumerate() {
                let replacement = if direct_structural_port(arg) {
                    if let Some(position) = matching_interface_position(arg, slots, &state.claimed)
                    {
                        state.claimed[position] = true;
                        if let Some(replacement) = replacements.get(&position) {
                            // Even an identical replacement must be recorded in
                            // every branch of a sum before unchanged atoms are reused.
                            state.applied[position] = true;
                            if replacement.as_view() == arg {
                                arg.into()
                            } else {
                                replacement.clone().into()
                            }
                        } else {
                            arg.into()
                        }
                    } else {
                        // Contracted dummy slots and representation-valued metadata
                        // are not part of the surviving public interface.
                        arg.into()
                    }
                } else if tensor_leaf || (fun.get_symbol() == SPENSO_TAG.chain && position < 2) {
                    // Tensor metadata and supplied chain endpoints are opaque
                    // to public port rewriting. Their compact representations
                    // do not denote surviving external axes.
                    arg.into()
                } else if is_composite_head(fun.get_symbol(), &SPENSO_TAG)
                    && !arg.is_tensorial(StrictTensorFilter::Tagged)
                    && !matches!(
                        arg,
                        AtomView::Fun(nested)
                            if nested.get_symbol() == *SYM
                                || nested.get_symbol() == *ANTISYM
                                || nested.get_symbol() == *CYCLIC
                    )
                {
                    arg.into()
                } else {
                    rewrite_ports(arg, slots, replacements, state)
                };
                rewritten.push(replacement);
            }
            if !value.needs_normalization()
                && rewritten
                    .iter()
                    .all(|arg| matches!(arg, AtomOrView::View(_)))
                && fun.get_symbol().get_evaluation_info().is_none()
                && (fun.get_symbol().get_normalization_function().is_none()
                    || (fun.get_symbol() == ETS.metric && fun.get_nargs() == 2 && {
                        // MetricNormalization can call an arbitrary vector normalizer
                        // only for a mixed explicit-slot / other-argument pair.
                        let mut matcher = SlotMatcher::default();
                        let mut args = fun.iter();
                        !matches!(
                            (
                                matcher.classify(args.next().unwrap()),
                                matcher.classify(args.next().unwrap())
                            ),
                            (SlotMatch::Explicit(_), SlotMatch::Other)
                                | (SlotMatch::Other, SlotMatch::Explicit(_))
                        )
                    }))
            {
                value.into()
            } else {
                // Rebuilding callback-bearing heads remains observable even
                // when their arguments, or the callback result, are unchanged.
                let result = FunctionBuilder::new(fun.get_symbol())
                    .add_args(rewritten)
                    .finish();
                result.into()
            }
        }
        _ => value.into(),
    }
}

pub(crate) fn rewrite_interface_ports(
    value: &SymbolicTensor<PartialStructure>,
    replacements: &HashMap<usize, Atom>,
    indices: Option<&ParseState<AbstractIndex>>,
) -> Result<Atom, TensorCompositionError> {
    if value.expression.as_view().is_zero() {
        return Ok(Atom::Zero);
    }

    let slots = value.structure.logical_slots();
    for &position in replacements.keys() {
        validate_position(&slots, position)?;
    }
    if value.borrowed_parts().product_is_normalized()
        && replacements.iter().all(|(&position, replacement)| {
            matching_interface_position(
                replacement.as_view(),
                &slots[position..=position],
                &[false],
            )
            .is_some()
        })
    {
        return Ok(value.expression.clone());
    }
    let mut state = PortRewriteState::new(slots.len());
    state.explicit_copies =
        value.proofs.validated.get() == Some(&true) && value.structure.open_positions().is_empty();
    if replacements.values().any(|replacement| {
        matches!(replacement.as_view(), AtomView::Fun(vector)
        if vector.get_symbol().has_tag(&SPENSO_TAG.rank1))
    }) {
        state.indices = Some(
            indices
                .cloned()
                .unwrap_or_else(|| SymbolicTensor::<PartialStructure>::reserved_dummies([value])),
        );
    }
    // Identity rewrites still walk every summand, but need no regrouping proof.
    state.bulk_arithmetic = replacements.iter().any(|(&position, replacement)| {
        matching_interface_position(replacement.as_view(), &slots[position..=position], &[false])
            .is_none()
    }) && PortRewriteState::exact_arithmetic(value.expression.as_view())
        && replacements
            .values()
            .all(|replacement| PortRewriteState::exact_arithmetic(replacement.as_view()));
    let atom =
        rewrite_ports(value.expression.as_view(), &slots, replacements, &mut state).into_owned();
    if !(value.proofs.validated.get() == Some(&true) && value.structure.open_positions().is_empty())
        && let Some(position) = replacements
            .keys()
            .copied()
            .find(|position| !state.applied[*position])
    {
        return Err(TensorCompositionError::MissingInterfacePort { position });
    }
    Ok(atom)
}

fn collect_interface_positions(
    value: AtomView<'_>,
    slots: &[PartialSlot],
    claimed: &mut [bool],
    positions: &mut Vec<usize>,
) {
    if direct_structural_port(value) {
        if let Some(position) = matching_interface_position(value, slots, claimed) {
            claimed[position] = true;
            positions.push(position);
        }
        return;
    }

    match value {
        // Every summand has the same public interface. Inspecting the first one
        // avoids collecting the same external ports repeatedly.
        AtomView::Add(add) => {
            if let Some(term) = add.iter().next() {
                collect_interface_positions(term, slots, claimed, positions);
            }
        }
        AtomView::Mul(mul) => {
            for factor in mul.iter() {
                let tensorial = factor.is_tensorial(StrictTensorFilter::Tagged)
                    || matches!(
                        factor,
                        AtomView::Fun(fun)
                            if fun.get_symbol() == *shadowing::SYM
                                || fun.get_symbol() == *shadowing::ANTISYM
                                || fun.get_symbol() == *shadowing::CYCLIC
                    );
                if tensorial {
                    collect_interface_positions(factor, slots, claimed, positions);
                }
            }
        }
        AtomView::Pow(pow) => {
            let (base, _) = pow.get_base_exp();
            collect_interface_positions(base, slots, claimed, positions);
        }
        // A compact dot has no public ports of its own.
        AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.dot => {}
        AtomView::Fun(fun) => {
            let symbol = fun.get_symbol();
            let tensor_leaf = is_tensor_leaf_head(symbol);
            let skip = usize::from(symbol == SPENSO_TAG.trace);
            for argument in fun.iter().skip(skip) {
                let tensor_argument = !tensor_leaf
                    && is_composite_head(symbol, &SPENSO_TAG)
                    && (argument.is_tensorial(StrictTensorFilter::Tagged)
                        || matches!(
                            argument,
                            AtomView::Fun(nested)
                                if nested.get_symbol() == *shadowing::SYM
                                    || nested.get_symbol() == *shadowing::ANTISYM
                                    || nested.get_symbol() == *shadowing::CYCLIC
                        ));
                if direct_structural_port(argument) || tensor_argument {
                    collect_interface_positions(argument, slots, claimed, positions);
                }
            }
        }
        _ => {}
    }
}

/// Order an established interface by its port occurrences in the source atom.
///
/// The caller must establish the interface and its occurrence correspondence.
/// This only reorders those ports; it does not infer or validate the tensor.
/// Unresolved ports retain the existing positional matching convention. Zero
/// and scalar interfaces retain their declared layout without inspecting syntax.
///
/// Cyclic trace normalization may rotate factors, so retaining the earlier
/// interface order would make `as_tensor(to_expression())` observe a different
/// logical interface.
///
/// See the [materialization example](SymbolicTensor::materialized).
///
/// # Errors
///
/// Returns [`TensorCompositionError::MissingInterfacePort`] if an established
/// port cannot be matched to its occurrence in the source expression.
pub fn interface_in_atom_order(
    atom: &Atom,
    interface: &PartialStructure,
) -> Result<PartialStructure, TensorCompositionError> {
    if atom.as_view().is_zero() || interface.canonical().is_scalar() {
        return Ok(interface.clone());
    }
    let slots = interface.logical_slots();
    let mut claimed = vec![false; slots.len()];
    let mut positions = Vec::with_capacity(slots.len());
    collect_interface_positions(atom.as_view(), &slots, &mut claimed, &mut positions);
    if let Some(position) = claimed.iter().position(|claimed| !claimed) {
        return Err(TensorCompositionError::MissingInterfacePort { position });
    }
    Ok(PartialStructure::from_logical_slots(
        positions.into_iter().map(|position| slots[position]),
    ))
}

fn without_positions(interface: &PartialStructure, positions: &[usize]) -> PartialStructure {
    PartialStructure::from_logical_slots(
        interface
            .logical_slots()
            .into_iter()
            .enumerate()
            .filter(|(position, _)| !positions.contains(position))
            .map(|(_, slot)| slot),
    )
}

fn concatenate_interfaces(
    left: impl IntoIterator<Item = PartialSlot>,
    right: impl IntoIterator<Item = PartialSlot>,
) -> PartialStructure {
    PartialStructure::from_logical_slots(left.into_iter().chain(right))
}

pub fn port_atom(slot: PartialSlot) -> Atom {
    match slot.aind {
        PartialIndex::Explicit(index) => slot.rep().slot::<AbstractIndex, _>(index).to_atom(),
        PartialIndex::Open(_) => slot.rep().to_symbolic([]),
    }
}

fn channel_factor(
    value: &SymbolicTensor<PartialStructure>,
    channel: MatrixChannel,
) -> Result<Atom, TensorCompositionError> {
    let replacements = HashMap::from([
        (channel.input, Atom::var(SPENSO_TAG.chain_in)),
        (channel.output, Atom::var(SPENSO_TAG.chain_out)),
    ]);
    rewrite_interface_ports(value, &replacements, None)
}

#[derive(Clone, Debug)]
pub enum ProductPlan {
    Scalar,
    Outer,
    Contract(Vec<PortPair>),
    Compose(MatrixChannel, MatrixChannel),
}

impl ProductPlan {
    /// Allocate all graph-port assignments for this already-selected product.
    /// Semantic interfaces choose pairings; materialized views use the same
    /// original logical positions and never infer pairings from fresh labels.
    #[allow(clippy::type_complexity)]
    pub fn alignment_indices(
        &self,
        left: &SymbolicTensor<PartialStructure>,
        right: &SymbolicTensor<PartialStructure>,
    ) -> Result<
        (HashMap<usize, AbstractIndex>, HashMap<usize, AbstractIndex>),
        TensorCompositionError,
    > {
        match self {
            Self::Scalar | Self::Outer => left.aligned_indices(right, &[]),
            Self::Contract(pairs) => left.aligned_indices(right, pairs),
            Self::Compose(left_channel, right_channel) => left.aligned_indices(
                right,
                &[PortPair {
                    left: left_channel.output,
                    right: right_channel.input,
                }],
            ),
        }
    }

    /// Apply a chosen positional composition to semantic or materialized ports.
    /// Relabeling an operand does not require choosing its pairings again.
    /// See the [example](SymbolicTensor::compose).
    ///
    /// # Errors
    ///
    /// Rejects operands that do not satisfy this plan, invalid ports, incompatible
    /// representations, ambiguous contractions, or invalid index multiplicities.
    pub fn apply(
        &self,
        left: &SymbolicTensor<PartialStructure>,
        right: &SymbolicTensor<PartialStructure>,
    ) -> Result<SymbolicTensor<PartialStructure>, TensorCompositionError> {
        match self {
            Self::Scalar => {
                let interface = if left.is_scalar() {
                    right.structure.clone()
                } else if right.is_scalar() {
                    left.structure.clone()
                } else {
                    return Err(TensorCompositionError::InvalidResultInterface(
                        "scalar multiplication requires a scalar operand".into(),
                    ));
                };
                let atom = left.expression.as_ref() * right.expression.as_ref();
                if !left.structure.open_positions().is_empty()
                    || !right.structure.open_positions().is_empty()
                {
                    validate_explicit_index_occurrences(&atom)?;
                }
                let value = if left.borrowed_parts().product_is_normalized()
                    && right.borrowed_parts().product_is_normalized()
                {
                    SymbolicTensor::from_normalized_parts(atom, interface)
                } else {
                    SymbolicTensor::new(atom, interface)
                };
                value.validate_rewrite("multiply", &[left, right])
            }
            Self::Outer => left.outer(right),
            Self::Contract(pairs) => left.contract_ports(right, pairs),
            Self::Compose(left_channel, right_channel) => {
                left.compose(right, *left_channel, *right_channel)
            }
        }
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Assign graph identities to unresolved ports and apply all explicit
    /// assignments simultaneously in the original logical interface.
    pub fn graph_port_indices(
        &self,
        replacements: &HashMap<usize, AbstractIndex>,
    ) -> Result<HashMap<usize, AbstractIndex>, TensorCompositionError> {
        let slots = self.structure.logical_slots();
        let owner = AbstractIndex::fresh_open_owner();
        let mut indices = self
            .structure
            .open_positions()
            .into_iter()
            .enumerate()
            .map(|(axis, position)| (position, AbstractIndex::Open { owner, axis }))
            .collect::<HashMap<_, _>>();
        for (&position, &index) in replacements {
            validate_position(&slots, position)?;
            indices.insert(position, index);
        }
        Ok(indices)
    }

    /// Plan every index assignment before either operand is rewritten.
    /// Reserve both complete sources and all existing interface labels once,
    /// including names which occur only in a later pair or opaque metadata.
    #[allow(clippy::type_complexity)]
    pub fn aligned_indices(
        &self,
        right: &Self,
        pairs: &[PortPair],
    ) -> Result<
        (HashMap<usize, AbstractIndex>, HashMap<usize, AbstractIndex>),
        TensorCompositionError,
    > {
        let left_slots = self.structure.logical_slots();
        let right_slots = right.structure.logical_slots();
        let targets = pairs
            .iter()
            .map(|&pair| {
                let left = validate_position(&left_slots, pair.left)?;
                let right = validate_position(&right_slots, pair.right)?;
                shared_index(left, right, pair).map(|index| (pair, index))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut left = self.graph_port_indices(&HashMap::new())?;
        let mut right_indices = right.graph_port_indices(&HashMap::new())?;
        let dummies = std::cell::OnceCell::new();
        for (pair, target) in targets {
            let target = target.unwrap_or_else(|| {
                dummies
                    .get_or_init(|| Self::reserved_dummies([self, right]))
                    .fresh_index()
            });
            left.insert(pair.left, target);
            right_indices.insert(pair.right, target);
        }
        Ok((left, right_indices))
    }

    /// Restrict semantic and materialized interfaces to the same surviving
    /// graph ports, preserving their shared original logical positions.
    /// See the [example](Self::materialized).
    ///
    /// # Errors
    ///
    /// Rejects unequal semantic/materialized ranks, unresolved materialized ports,
    /// external graph ports missing from the materialized interface, or invalid chain closure.
    pub fn retain_materialized_ports(
        self,
        materialized: Self,
        external: &[LibrarySlot<AbstractIndex>],
    ) -> Result<(Self, Self), TensorCompositionError> {
        let semantic_slots = self.structure.logical_slots();
        let materialized_slots = materialized.structure.logical_slots();
        if semantic_slots.len() != materialized_slots.len() {
            return Err(TensorCompositionError::InvalidResultInterface(
                "semantic and materialized interfaces have different ranks".into(),
            ));
        }
        let mut dangling = external.to_vec();
        let mut retained = Vec::new();
        for (position, slot) in materialized_slots.iter().enumerate() {
            let PartialIndex::Explicit(index) = slot.aind else {
                return Err(TensorCompositionError::InvalidResultInterface(
                    "materialized interface contains an unresolved port".into(),
                ));
            };
            let slot = slot.rep().slot(index).to_lib();
            if let Some(found) = dangling.iter().position(|candidate| *candidate == slot) {
                dangling.swap_remove(found);
                retained.push(position);
            }
        }
        if !dangling.is_empty() {
            return Err(TensorCompositionError::InvalidResultInterface(format!(
                "graph ports are absent from the materialized interface: {dangling:?}"
            )));
        }
        let mut semantic = Self::from_normalized_parts(
            self.expression,
            PartialStructure::from_logical_slots(
                retained.iter().map(|&position| semantic_slots[position]),
            ),
        );
        semantic.proofs = self.proofs.algebra_result();
        let semantic = semantic.normalize_closed_root_chain()?;
        let mut retained_materialized = Self::from_normalized_parts(
            materialized.expression,
            PartialStructure::from_logical_slots(
                retained
                    .iter()
                    .map(|&position| materialized_slots[position]),
            ),
        );
        retained_materialized.proofs = materialized.proofs.algebra_result();
        let retained_materialized = retained_materialized.normalize_closed_root_chain()?;
        Ok((semantic, retained_materialized))
    }

    pub fn trace_indices(
        &self,
        channel: MatrixChannel,
    ) -> Result<HashMap<usize, AbstractIndex>, TensorCompositionError> {
        let pair = PortPair {
            left: channel.input,
            right: channel.output,
        };
        let target = self.contraction_index(self, pair)?;
        self.graph_port_indices(&HashMap::from([(pair.left, target), (pair.right, target)]))
    }

    fn validate_rewrite(
        self,
        _operation: &str,
        sources: &[&Self],
    ) -> Result<Self, TensorCompositionError> {
        if sources
            .iter()
            .all(|source| source.structure.open_positions().is_empty())
        {
            if let Some(occurrences) = self.proofs.occurrences.get() {
                occurrences.validate()?;
            } else {
                let occurrences = ExplicitIndexOccurrences::product(sources)?;
                let _ = self.proofs.occurrences.set(Arc::new(occurrences));
            }
            #[cfg(debug_assertions)]
            TensorProofs::assert_normalization_preserves_interface(
                _operation,
                &sources
                    .iter()
                    .map(|source| &source.expression)
                    .collect::<Vec<_>>(),
                &self.expression,
                Ok(&self.structure),
                None,
            );
            let _ = self.proofs.validated.set(true);
            if sources
                .iter()
                .all(|source| source.proofs.intrinsic.get() == Some(&true))
            {
                let _ = self.proofs.intrinsic.set(true);
            }
            if sources
                .iter()
                .all(|source| source.borrowed_parts().product_is_normalized())
            {
                let _ = self.proofs.normalized_product.set(true);
            }
            return Ok(self);
        }
        if sources
            .iter()
            .any(|source| !source.normalization_is_intrinsic())
        {
            // Materialization can make a pair explicit before its enclosing
            // composition removes it. Compare the contracted external interface.
            let external = InterfaceInference::merge_explicit_interface_sequence(slice::from_ref(
                &self.structure,
            ))
            .map_err(TensorCompositionError::from)?;
            Self::validate_interface(&self.expression, &external)
                .map_err(TensorCompositionError::from)?;
        }
        let mut occurrences = match self.proofs.occurrences.get() {
            Some(occurrences) => occurrences.as_ref().clone(),
            None => ExplicitIndexOccurrences::product(sources)?,
        };
        occurrences.branch_max(&ExplicitIndexOccurrences::from_interface(&self.structure));
        occurrences.validate()?;
        let _ = self.proofs.occurrences.set(Arc::new(occurrences));
        if sources
            .iter()
            .all(|source| source.proofs.validated.get() == Some(&true))
        {
            let _ = self.proofs.validated.set(true);
        }
        Ok(self)
    }

    /// See the [example](Self::materialized).
    ///
    /// # Errors
    ///
    /// Rejects axes that are not a rank-preserving permutation or invalid source metadata.
    pub fn permuted(&self, axes: &[usize]) -> Result<Self, TensorCompositionError> {
        let rank = self.rank();
        let mut seen = vec![false; rank];
        if axes.len() != rank
            || axes
                .iter()
                .any(|&axis| axis >= rank || mem::replace(&mut seen[axis], true))
        {
            return Err(TensorCompositionError::InvalidAxisPermutation {
                rank,
                axes: axes.to_vec(),
            });
        }
        self.explicit_index_occurrences()
            .map_err(TensorCompositionError::from)?;
        if self.structure.open_positions().is_empty() {
            if axes.iter().copied().eq(0..rank) {
                return Ok(self.clone());
            }
            let slots = self.structure.logical_slots();
            let mut value = self.clone();
            value.structure =
                PartialStructure::from_logical_slots(axes.iter().map(|&axis| slots[axis]));
            value.proofs = self.proofs.algebra_result();
            return Ok(value);
        }
        // Give unresolved occurrences identities before changing positional order.
        let dummies = OnceCell::new();
        let replacements = self
            .structure
            .open_positions()
            .into_iter()
            .map(|position| {
                (
                    position,
                    dummies
                        .get_or_init(|| Self::reserved_dummies([self]))
                        .fresh_index(),
                )
            })
            .collect();
        let mut value = self.reindex_interface_ports(&replacements)?;
        let slots = value.structure.logical_slots();
        value.proofs = value.proofs.algebra_result();
        value.structure =
            PartialStructure::from_logical_slots(axes.iter().map(|&axis| slots[axis]));
        Ok(value)
    }

    /// Give each unresolved public port a fresh explicit identity, preserving logical order.
    ///
    /// # Errors
    ///
    /// Returns an error if admission, placeholder scopes, or index multiplicities are invalid.
    ///
    /// # Example
    ///
    /// ```
    /// use std::collections::HashMap;
    /// use std::error::Error;
    /// use idenso::tensor::{SymbolicTensor, composition};
    /// use idenso::tensor::composition::PortPair;
    /// use spenso::structure::abstract_index::AbstractIndex;
    /// use spenso::structure::representation::{Euclidean, RepName};
    /// use spenso::tensor_symbol;
    /// use symbolica::function;
    ///
    /// let rep = Euclidean {}.new_rep(3).to_lib();
    /// let open = SymbolicTensor::infer(function!(tensor_symbol!("docs::open"), rep.to_symbolic([])))?;
    /// assert_eq!(open.materialized()?.rank(), 1);
    /// let index = AbstractIndex::Normal(11);
    /// let assignments = HashMap::from([(0, index)]);
    /// let expression = open.materialize_interface_ports(&assignments)?;
    /// let stored = open.reindex_interface_ports(&assignments)?;
    /// assert_eq!(stored.expression(), &expression);
    /// assert!(stored.permuted(&[0])?.same_value(&stored));
    /// let layout = composition::interface_in_atom_order(stored.expression(), stored.structure())?;
    /// assert_eq!(&layout, stored.structure());
    /// assert_eq!(open.outer(&open)?.rank(), 2);
    /// assert!(open.contract_ports(&open, &[PortPair { left: 0, right: 0 }])?.is_scalar());
    /// let (semantic, stored) = open.retain_materialized_ports(stored, &[rep.slot(index)])?;
    /// assert_eq!(semantic.rank(), stored.rank());
    /// # Ok::<(), Box<dyn Error>>(())
    /// ```
    pub fn materialized(&self) -> Result<Self, TensorCompositionError> {
        if self.structure.open_positions().is_empty() {
            self.explicit_index_occurrences()
                .map_err(TensorCompositionError::from)?;
            return Ok(self.clone());
        }
        let dummies = OnceCell::new();
        let replacements = self
            .structure
            .open_positions()
            .into_iter()
            .map(|position| {
                (
                    position,
                    dummies
                        .get_or_init(|| Self::reserved_dummies([self]))
                        .fresh_index(),
                )
            })
            .collect::<HashMap<_, _>>();
        self.reindex_interface_ports(&replacements)
    }

    pub fn matrix_channel(&self) -> Option<MatrixChannel> {
        let slots = self.structure.logical_slots();
        if slots.len() >= 2
            && has_transparent_chain_channel(self.expression.as_view(), slots[0], slots[1])
        {
            return Some(MatrixChannel {
                input: 0,
                output: 1,
            });
        }

        let pairs = (0..slots.len())
            .flat_map(|left| ((left + 1)..slots.len()).map(move |right| (left, right)))
            .filter(|&(left, right)| representations_match(&slots[left], &slots[right]))
            .collect::<Vec<_>>();
        let [(first, second)] = pairs.as_slice() else {
            return None;
        };

        let first_rep = representation(&slots[*first]);
        let second_rep = representation(&slots[*second]);
        if first_rep.rep.is_self_dual() || (first_rep.rep.is_base() && second_rep.rep.is_dual()) {
            Some(MatrixChannel {
                input: *first,
                output: *second,
            })
        } else {
            Some(MatrixChannel {
                input: *second,
                output: *first,
            })
        }
    }

    /// Replace selected unresolved ports by explicit indices in every summand.
    ///
    /// Positions are interpreted in the public logical interface, never in the
    /// canonical storage order of `OrderedStructure`.
    /// See the [example](Self::materialized).
    ///
    /// # Errors
    ///
    /// Rejects invalid positions, changing an existing explicit index, missing source
    /// ports, or invalid resulting index multiplicities.
    pub fn materialize_interface_ports(
        &self,
        replacements: &HashMap<usize, AbstractIndex>,
    ) -> Result<Atom, TensorCompositionError> {
        let slots = self.structure.logical_slots();
        if self.borrowed_parts().product_is_normalized()
            && replacements.iter().all(|(&position, &index)| {
                slots
                    .get(position)
                    .is_some_and(|slot| slot.aind == PartialIndex::Explicit(index))
            })
        {
            return Ok(self.expression.clone());
        }
        let replacements = replacements
            .iter()
            .map(|(&position, &index)| {
                let slot = validate_position(&slots, position)?;
                let atom = match slot.aind {
                    PartialIndex::Explicit(current) if current == index => port_atom(slot),
                    PartialIndex::Explicit(_) => {
                        return Err(TensorCompositionError::UnequalExplicitIndices {
                            left: position,
                            right: position,
                        });
                    }
                    PartialIndex::Open(_) => slot.rep().slot::<AbstractIndex, _>(index).to_atom(),
                };
                Ok((position, atom))
            })
            .collect::<Result<HashMap<_, _>, TensorCompositionError>>()?;
        let atom = rewrite_interface_ports(self, &replacements, None)?;
        validate_explicit_index_occurrences(&atom)?;
        Ok(atom)
    }

    /// Replace selected public ports by explicit graph indices.
    ///
    /// Unlike `materialize_interface_ports`, this also permits changing an
    /// already-explicit index. It is used to keep a tensor network's execution
    /// graph collision-free while its semantic interface remains unchanged.
    /// See the [example](Self::materialized).
    ///
    /// # Errors
    ///
    /// Rejects invalid positions, missing source ports, or a relabeling that introduces
    /// more than two compatible uses of an index. Debug builds panic if a user
    /// normalizer violates its free-index preservation contract.
    pub fn reindex_interface_ports(
        &self,
        replacements: &HashMap<usize, AbstractIndex>,
    ) -> Result<SymbolicTensor<PartialStructure>, TensorCompositionError> {
        let explicit = self.structure.open_positions().is_empty();
        let mut occurrences = ExplicitIndexOccurrences::product(&[self])?;
        let mut logical = self.structure.logical_slots();
        let mut mapped = Vec::new();
        let mut materialized = Vec::new();
        let replacements = replacements
            .iter()
            .map(|(&position, &index)| {
                let slot = validate_position(&logical, position)?;
                logical[position].set_aind(PartialIndex::Explicit(index));
                if let PartialIndex::Explicit(previous) = slot.aind {
                    mapped.push((slot.rep().slot(previous), slot.rep().slot(index)));
                } else {
                    materialized.push(slot.rep().slot(index));
                }
                Ok((
                    position,
                    slot.rep().slot::<AbstractIndex, _>(index).to_atom(),
                ))
            })
            .collect::<Result<HashMap<_, _>, TensorCompositionError>>()?;
        if explicit
            && self.borrowed_parts().product_is_normalized()
            && mapped.iter().all(|(source, target)| source == target)
        {
            return Ok(self.clone());
        }
        let atom = rewrite_interface_ports(self, &replacements, None)?;
        // Decrement every old use first, so swaps and partial relabellings of
        // repeated free indices do not create transient third occurrences.
        for &(source, _) in &mapped {
            occurrences.remove(source)?;
        }
        for &(_, target) in &mapped {
            occurrences.add(target, 1);
        }
        for slot in materialized {
            occurrences.add(slot, 1);
        }
        // If relabelling splits a repeated key, the surviving public singleton
        // supplies its exact orientation rather than the old pair's representative.
        for slot in &logical {
            if let PartialIndex::Explicit(index) = slot.aind {
                let slot = slot.rep().slot(index);
                if occurrences.count(slot) == 1 {
                    occurrences.remove(slot)?;
                    occurrences.add(slot, 1);
                }
            }
        }
        occurrences.validate()?;
        let interface = PartialStructure::from_logical_slots(logical);
        let mut value = if explicit {
            Self::from_normalized_parts(atom, interface)
        } else {
            validate_explicit_index_occurrences(&atom)?;
            Self::new(atom, interface)
        };
        if explicit {
            value.proofs = self.proofs.algebra_result();
            if mapped
                .iter()
                .any(|(_, slot)| matches!(slot.aind(), AbstractIndex::Open { .. }))
            {
                value.proofs.normalized_product = Default::default();
            }
        }
        value.proofs.occurrences = OnceLock::from(Arc::new(occurrences));
        if explicit {
            TensorProofs::assert_normalization_preserves_interface(
                "reindex_interface_ports",
                &[&self.expression],
                &value.expression,
                Ok(&value.structure),
                None,
            );
            Ok(value)
        } else {
            value.validate_rewrite("reindex_interface_ports", &[self])
        }
    }

    /// Canonicalize a root chain whose two explicit endpoints have become one contraction.
    /// See the [example](Self::compose).
    ///
    /// # Errors
    ///
    /// Rejects an inconsistent chain interface, missing source ports, or invalid
    /// explicit-index multiplicities while closing the root channel.
    pub fn normalize_closed_root_chain(
        self,
    ) -> Result<SymbolicTensor<PartialStructure>, TensorCompositionError> {
        let Some((coefficient, function)) = root_chain_parts(self.expression.as_view()) else {
            return Ok(self);
        };
        let arguments = function.iter().collect::<Vec<_>>();
        let [start, end, factors @ ..] = arguments.as_slice() else {
            return Ok(self);
        };
        let (Ok(start), Ok(end)) = (
            Slot::<LibraryRep, AbstractIndex>::try_from(*start),
            Slot::<LibraryRep, AbstractIndex>::try_from(*end),
        ) else {
            return Ok(self);
        };
        let start_rep = start.rep();
        let end_rep = end.rep();
        if start.aind() != end.aind()
            || !start_rep.matches(&end_rep)
            || !(start_rep.rep.is_self_dual() || start_rep.rep.is_base() && end_rep.rep.is_dual())
        {
            return Ok(self);
        }

        let channel = MatrixChannel {
            input: 0,
            output: 1,
        };
        let interface = if root_chain_channel_is_live(&self, channel) {
            without_positions(&self.structure, &[channel.input, channel.output])
        } else {
            self.structure.clone()
        };
        let atom = coefficient
            * shadowing::trace(
                start_rep.to_symbolic([]),
                factors.iter().map(|factor| factor.to_owned()),
            );
        let interface = if self.structure.open_positions().is_empty() {
            interface
        } else {
            interface_in_atom_order(&atom, &interface)?
        };
        SymbolicTensor::from_normalized_parts(atom, interface)
            .validate_rewrite("normalize_closed_root_chain", &[&self])
    }

    /// See the [example](Self::materialized).
    ///
    /// # Errors
    ///
    /// Rejects equal explicit indices on compatible operand ports, invalid placeholder
    /// scopes, or invalid index multiplicities.
    pub fn outer(
        &self,
        right: &SymbolicTensor<PartialStructure>,
    ) -> Result<Self, TensorCompositionError> {
        let left_slots = self.structure.logical_slots();
        let right_slots = right.structure.logical_slots();
        let collisions = compatible_pairs(&self.structure, &right.structure)
            .into_iter()
            .filter(|pair| {
                matches!(left_slots[pair.left].aind, PartialIndex::Explicit(_))
                    && matches!(right_slots[pair.right].aind, PartialIndex::Explicit(_))
            })
            .collect::<Vec<_>>();
        if !collisions.is_empty() {
            return Err(TensorCompositionError::ExplicitOuterContraction(collisions));
        }
        let atom = if self.expression.as_view().is_zero() || right.expression.as_view().is_zero() {
            Atom::Zero
        } else {
            SymbolicTensor::bracket_product(self.expression.as_view(), right.expression.as_view())
        };
        if !self.structure.open_positions().is_empty()
            || !right.structure.open_positions().is_empty()
        {
            validate_explicit_index_occurrences(&atom)?;
        }
        SymbolicTensor::from_normalized_parts(
            atom,
            concatenate_interfaces(
                self.structure.logical_slots(),
                right.structure.logical_slots(),
            ),
        )
        .validate_rewrite("outer", &[self, right])
    }

    /// See the [example](Self::materialized).
    ///
    /// # Errors
    ///
    /// Rejects invalid or incompatible port pairs, unequal explicit indices, missing
    /// source ports, unsupported nested chains, or excessive index multiplicity.
    pub fn contract_ports(
        &self,
        right: &SymbolicTensor<PartialStructure>,
        positions: &[PortPair],
    ) -> Result<SymbolicTensor<PartialStructure>, TensorCompositionError> {
        if positions.is_empty() {
            return self.outer(right);
        }
        let left_slots = self.structure.logical_slots();
        let right_slots = right.structure.logical_slots();
        let mut pairs = positions.to_vec();
        pairs.extend(
            left_slots
                .iter()
                .enumerate()
                .filter_map(|(left_position, left_slot)| {
                    (!positions.iter().any(|pair| pair.left == left_position))
                        .then_some((left_position, left_slot))
                })
                .flat_map(|(left_position, left_slot)| {
                    right_slots.iter().enumerate().filter_map(
                        move |(right_position, right_slot)| {
                            (!positions.iter().any(|pair| pair.right == right_position)
                            && matches!(
                                (left_slot.aind, right_slot.aind),
                                (PartialIndex::Explicit(left), PartialIndex::Explicit(right))
                                    if left == right && representations_match(left_slot, right_slot)
                            ))
                        .then_some(PortPair {
                            left: left_position,
                            right: right_position,
                        })
                        },
                    )
                }),
        );
        let mut contracted_left = HashSet::new();
        let mut contracted_right = HashSet::new();
        for pair in &pairs {
            if !contracted_left.insert(pair.left) || !contracted_right.insert(pair.right) {
                return Err(TensorCompositionError::Ambiguous(pairs));
            }
        }
        let indices = pairs
            .iter()
            .map(|&pair| {
                let left_slot = validate_position(&left_slots, pair.left)?;
                let right_slot = validate_position(&right_slots, pair.right)?;
                shared_index(left_slot, right_slot, pair).map(|index| (pair, index))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let interface = PartialStructure::from_logical_slots(
            left_slots
                .iter()
                .copied()
                .enumerate()
                .filter(|(position, _)| !contracted_left.contains(position))
                .map(|(_, slot)| slot)
                .chain(
                    right_slots
                        .iter()
                        .copied()
                        .enumerate()
                        .filter(|(position, _)| !contracted_right.contains(position))
                        .map(|(_, slot)| slot),
                ),
        );
        let dummies = OnceCell::new();
        let mut left_indices = HashMap::new();
        let mut right_indices = HashMap::new();
        for (pair, shared) in indices {
            let index = shared.unwrap_or_else(|| {
                dummies
                    .get_or_init(|| Self::reserved_dummies([self, right]))
                    .fresh_index()
            });
            left_indices.insert(pair.left, index);
            right_indices.insert(pair.right, index);
        }
        let explicit = self.structure.open_positions().is_empty()
            && right.structure.open_positions().is_empty();
        let materialized;
        let (materialized_left, materialized_right) = if explicit {
            (self, right)
        } else {
            materialized = (
                self.reindex_interface_ports(&left_indices)?,
                right.reindex_interface_ports(&right_indices)?,
            );
            (&materialized.0, &materialized.1)
        };
        let atom = if self.is_zero() || right.is_zero() {
            Atom::Zero
        } else if self.rank() == 1
            && right.rank() == 1
            && [self, right].into_iter().all(|operand| {
                matches!(operand.expression.as_view(), AtomView::Fun(vector)
                    if vector.get_symbol().has_tag(&SPENSO_TAG.rank1))
            })
        {
            // Only declared vector heads choose dot notation automatically.
            // A composite with one surviving axis retains its indexed factors;
            // the structural contractor can pack them into Schoonschip notation.
            self.compact_dot_expression(right)?
        } else {
            SymbolicTensor::bracket_product(
                materialized_left.expression.as_view(),
                materialized_right.expression.as_view(),
            )
        };
        if !explicit {
            validate_explicit_index_occurrences(&atom)?;
        }
        SymbolicTensor::from_normalized_parts(atom, interface)
            .validate_rewrite("contract_ports", &[materialized_left, materialized_right])
    }

    fn compact_dot_expression(&self, right: &Self) -> Result<Atom, TensorCompositionError> {
        let compact = |operand: &Self| {
            let slot = operand.structure.logical_slots()[0];
            rewrite_interface_ports(
                operand,
                &HashMap::from([(0, slot.rep().to_symbolic([]))]),
                None,
            )
        };
        Ok(FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(compact(self)?)
            .add_arg(compact(right)?)
            .finish())
    }

    pub fn chain_factors(
        &self,
        channel: MatrixChannel,
    ) -> Result<Vec<Atom>, TensorCompositionError> {
        let (coefficient, mut factors) = self.chain_parts(channel)?;
        if coefficient != Atom::num(1) {
            if let Some(first) = factors.first_mut() {
                *first *= coefficient;
            } else {
                factors.push(coefficient);
            }
        }
        Ok(factors)
    }

    fn chain_parts(
        &self,
        channel: MatrixChannel,
    ) -> Result<(Atom, Vec<Atom>), TensorCompositionError> {
        if self.structure.open_positions().is_empty() {
            self.explicit_index_occurrences()
                .map_err(TensorCompositionError::from)?;
        } else {
            // Raw positional interfaces can deliberately differ from a
            // builtin's registered argument order. The port rewrite below
            // owns their correspondence; only placeholder scopes are needed here.
            self.expression.validate_chain_like_nesting()?;
        }
        if root_chain_channel_is_live(self, channel) {
            let Some((coefficient, function)) = root_chain_parts(self.expression.as_view()) else {
                unreachable!("a live root chain channel requires a chain function")
            };
            let factors = function
                .iter()
                .skip(2)
                .map(|factor| factor.to_owned())
                .collect::<Vec<_>>();
            return Ok((coefficient, factors));
        }
        if contains_chain_like(self.expression.as_view()) {
            return Err(TensorCompositionError::NestedChainLike);
        }
        Ok((Atom::num(1), vec![channel_factor(self, channel)?]))
    }

    /// Join the selected output and input channels, retaining uncontracted spectator ports.
    ///
    /// # Errors
    ///
    /// Rejects invalid or incompatible channels, ambiguous spectator contractions,
    /// invalid index multiplicities, and unsupported nested chain scopes.
    ///
    /// # Example
    ///
    /// ```
    /// use std::error::Error;
    /// use idenso::tensor::SymbolicTensor;
    /// use idenso::tensor::composition::{MatrixChannel, ProductPlan};
    /// use spenso::structure::representation::{Euclidean, RepName};
    /// use spenso::tensor_symbol;
    /// use symbolica::function;
    ///
    /// let rep = Euclidean {}.new_rep(3).to_lib();
    /// let matrix = SymbolicTensor::infer(function!(
    ///     tensor_symbol!("docs::matrix"), rep.to_symbolic([]), rep.to_symbolic([])
    /// ))?;
    /// let channel = MatrixChannel { input: 0, output: 1 };
    /// let product = matrix.compose(&matrix, channel, channel)?;
    /// assert_eq!(product.rank(), 2);
    /// assert!(product.trace_ports(channel)?.is_scalar());
    /// assert!(product.trace_unique()?.is_scalar());
    /// assert_eq!(product.clone().normalize_closed_root_chain()?.rank(), 2);
    /// let planned = ProductPlan::Compose(channel, channel).apply(&matrix, &matrix)?;
    /// assert!(planned.same_value(&product));
    /// # Ok::<(), Box<dyn Error>>(())
    /// ```
    pub fn compose(
        &self,
        right: &SymbolicTensor<PartialStructure>,
        left_channel: MatrixChannel,
        right_channel: MatrixChannel,
    ) -> Result<SymbolicTensor<PartialStructure>, TensorCompositionError> {
        if left_channel.input == left_channel.output {
            return Err(TensorCompositionError::DegenerateChannel {
                input: left_channel.input,
                output: left_channel.output,
            });
        }
        if right_channel.input == right_channel.output {
            return Err(TensorCompositionError::DegenerateChannel {
                input: right_channel.input,
                output: right_channel.output,
            });
        }
        let left_slots = self.structure.logical_slots();
        let right_slots = right.structure.logical_slots();
        let left_input = validate_position(&left_slots, left_channel.input)?;
        let left_output = validate_position(&left_slots, left_channel.output)?;
        let right_input = validate_position(&right_slots, right_channel.input)?;
        let right_output = validate_position(&right_slots, right_channel.output)?;
        if !representations_match(&left_input, &left_output) {
            return Err(TensorCompositionError::IncompatiblePorts {
                left: left_channel.input,
                right: left_channel.output,
            });
        }
        let left_input_rep = representation(&left_input);
        let left_output_rep = representation(&left_output);
        if !(left_input_rep.rep.is_self_dual()
            || left_input_rep.rep.is_base() && left_output_rep.rep.is_dual())
        {
            return Err(TensorCompositionError::InvalidChannelOrientation {
                input: left_channel.input,
                output: left_channel.output,
            });
        }
        if !representations_match(&right_input, &right_output) {
            return Err(TensorCompositionError::IncompatiblePorts {
                left: right_channel.input,
                right: right_channel.output,
            });
        }
        let right_input_rep = representation(&right_input);
        let right_output_rep = representation(&right_output);
        if !(right_input_rep.rep.is_self_dual()
            || right_input_rep.rep.is_base() && right_output_rep.rep.is_dual())
        {
            return Err(TensorCompositionError::InvalidChannelOrientation {
                input: right_channel.input,
                output: right_channel.output,
            });
        }
        let shared = shared_index(
            left_output,
            right_input,
            PortPair {
                left: left_channel.output,
                right: right_channel.input,
            },
        )?;
        let spectator_pairs = left_slots
            .iter()
            .enumerate()
            .filter(|(position, _)| {
                *position != left_channel.input && *position != left_channel.output
            })
            .flat_map(|(left_position, left_slot)| {
                right_slots
                    .iter()
                    .enumerate()
                    .filter(move |(position, _)| {
                        *position != right_channel.input && *position != right_channel.output
                    })
                    .filter_map(move |(right_position, right_slot)| {
                        matches!(
                            (left_slot.aind, right_slot.aind),
                            (PartialIndex::Explicit(left), PartialIndex::Explicit(right))
                                if left == right && representations_match(left_slot, right_slot)
                        )
                        .then_some(PortPair {
                            left: left_position,
                            right: right_position,
                        })
                    })
            })
            .collect::<Vec<_>>();
        if spectator_pairs.iter().enumerate().any(|(position, pair)| {
            spectator_pairs[position + 1..]
                .iter()
                .any(|candidate| pair.left == candidate.left || pair.right == candidate.right)
        }) {
            return Err(TensorCompositionError::Ambiguous(spectator_pairs));
        }
        let mut occurrences = ExplicitIndexOccurrences::product(&[self, right])?;
        if let Some(index) = shared {
            occurrences.add(
                left_output.rep().slot(index),
                usize::from(matches!(left_output.aind, PartialIndex::Open(_)))
                    + usize::from(matches!(right_input.aind, PartialIndex::Open(_))),
            );
            occurrences.validate()?;
        }
        let contracted_left = spectator_pairs
            .iter()
            .map(|pair| pair.left)
            .collect::<HashSet<_>>();
        let contracted_right = spectator_pairs
            .iter()
            .map(|pair| pair.right)
            .collect::<HashSet<_>>();

        let interface = PartialStructure::from_logical_slots(
            [left_input, right_output]
                .into_iter()
                .chain(
                    left_slots
                        .into_iter()
                        .enumerate()
                        .filter(|(position, _)| {
                            *position != left_channel.input
                                && *position != left_channel.output
                                && !contracted_left.contains(position)
                        })
                        .map(|(_, slot)| slot),
                )
                .chain(
                    right_slots
                        .into_iter()
                        .enumerate()
                        .filter(|(position, _)| {
                            *position != right_channel.input
                                && *position != right_channel.output
                                && !contracted_right.contains(position)
                        })
                        .map(|(_, slot)| slot),
                ),
        );
        if self.expression.as_view().is_zero() || right.expression.as_view().is_zero() {
            let value = SymbolicTensor::from_normalized_parts(Atom::Zero, interface);
            let _ = value.proofs.occurrences.set(Arc::new(occurrences));
            return value.validate_rewrite("compose", &[self, right]);
        }

        let (coefficient, factors) = self.chain_parts(left_channel)?;
        let (right_coefficient, right_factors) = right.chain_parts(right_channel)?;
        let atom = coefficient
            * right_coefficient
            * SPENSO_TAG.chain(
                port_atom(left_input),
                port_atom(right_output),
                factors.into_iter().chain(right_factors),
            );
        if self.structure.open_positions().is_empty() && right.structure.open_positions().is_empty()
        {
            let value = SymbolicTensor::from_normalized_parts(atom, interface);
            let _ = value.proofs.occurrences.set(Arc::new(occurrences));
            value
                .validate_rewrite("compose", &[self, right])?
                .normalize_closed_root_chain()
        } else {
            validate_explicit_index_occurrences(&atom)?;
            let value = SymbolicTensor::new(atom, interface);
            let _ = value.proofs.occurrences.set(Arc::new(occurrences));
            value
                .normalize_closed_root_chain()?
                .validate_rewrite("compose", &[self, right])
        }
    }

    /// Choose explicit contractions first, then established matrix composition.
    /// Other pairings maximize the number of contractions and use the number of
    /// unresolved-unresolved pairs to break ties. Equally preferred plans are ambiguous.
    pub fn product_plan(
        &self,
        right: &SymbolicTensor<PartialStructure>,
    ) -> Result<ProductPlan, TensorCompositionError> {
        if self.is_scalar() || right.is_scalar() {
            return Ok(ProductPlan::Scalar);
        }
        let left_slots = self.structure.logical_slots();
        let right_slots = right.structure.logical_slots();
        let candidates = compatible_pairs(&self.structure, &right.structure);
        // Repeated explicit labels fix their partners before unresolved ports are considered.
        let mut pairs = candidates
            .iter()
            .copied()
            .filter(|pair| {
                matches!(
                    (left_slots[pair.left].aind, right_slots[pair.right].aind),
                    (PartialIndex::Explicit(_), PartialIndex::Explicit(_))
                )
            })
            .collect::<Vec<_>>();
        let mut named_left = HashSet::new();
        let mut named_right = HashSet::new();
        for pair in &pairs {
            if !named_left.insert(pair.left) || !named_right.insert(pair.right) {
                return Err(TensorCompositionError::Ambiguous(pairs));
            }
        }
        let remaining = candidates
            .iter()
            .copied()
            .filter(|pair| !named_left.contains(&pair.left) && !named_right.contains(&pair.right))
            .collect::<Vec<_>>();
        // Fully determined Einstein contractions need no matrix-channel inference.
        // A rank-two composite may expose its ports on separate vector factors,
        // which cannot each supply both placeholders required by a chain factor.
        if !remaining.is_empty()
            && let (Some(left_channel), Some(right_channel)) =
                (self.matrix_channel(), right.matrix_channel())
        {
            let channel_pair = PortPair {
                left: left_channel.output,
                right: right_channel.input,
            };
            // Matrix composition keeps its channel convention, but cannot override named contractions.
            if candidates.contains(&channel_pair)
                && pairs.iter().all(|pair| {
                    *pair == channel_pair
                        || (pair.left == left_channel.input && pair.right == right_channel.output)
                        || (pair.left != left_channel.input
                            && pair.left != left_channel.output
                            && pair.right != right_channel.input
                            && pair.right != right_channel.output)
                })
            {
                return Ok(ProductPlan::Compose(left_channel, right_channel));
            }
        }
        // Minimum-cost augmenting paths maximize simultaneous contractions first,
        // then prefer unresolved-unresolved pairs. Matched edges run backwards so
        // a later augmentation can revise an earlier choice without losing that priority.
        let pair_cost = |pair: PortPair| {
            isize::from(!matches!(
                (left_slots[pair.left].aind, right_slots[pair.right].aind),
                (PartialIndex::Open(_), PartialIndex::Open(_))
            ))
        };
        let best_matching = |excluded: Option<PortPair>| {
            let mut left_matches: Vec<Option<usize>> = vec![None; left_slots.len()];
            let mut right_matches: Vec<Option<usize>> = vec![None; right_slots.len()];
            let right_offset = left_slots.len();
            loop {
                let mut predecessors = vec![None; right_slots.len()];
                let mut distances = vec![None; left_slots.len() + right_slots.len()];
                for (position, matched) in left_matches.iter().enumerate() {
                    if matched.is_none() {
                        distances[position] = Some(0isize);
                    }
                }
                // Reverse edges can have negative costs, so use Bellman-Ford rather
                // than a breadth-first search. Optimal residual matchings have no negative cycles.
                for _ in 0..distances.len() {
                    let mut changed = false;
                    for &pair in remaining.iter().filter(|&&pair| Some(pair) != excluded) {
                        let (from, to, cost) = if left_matches[pair.left] == Some(pair.right) {
                            (right_offset + pair.right, pair.left, -pair_cost(pair))
                        } else {
                            (pair.left, right_offset + pair.right, pair_cost(pair))
                        };
                        let Some(distance) = distances[from] else {
                            continue;
                        };
                        let distance = distance + cost;
                        if distances[to].is_none_or(|current| distance < current) {
                            distances[to] = Some(distance);
                            if to >= right_offset {
                                predecessors[pair.right] = Some(pair.left);
                            }
                            changed = true;
                        }
                    }
                    if !changed {
                        break;
                    }
                }
                let target = right_matches
                    .iter()
                    .enumerate()
                    .filter(|(_, matched)| matched.is_none())
                    .filter_map(|(position, _)| {
                        distances[right_offset + position].map(|cost| (position, cost))
                    })
                    .min_by_key(|&(_, cost)| cost);
                let Some((mut target, _)) = target else {
                    break;
                };
                loop {
                    let source = predecessors[target].unwrap();
                    let previous = left_matches[source].replace(target);
                    right_matches[target] = Some(source);
                    match previous {
                        Some(previous) => target = previous,
                        None => break,
                    }
                }
            }
            right_matches
                .into_iter()
                .enumerate()
                .filter_map(|(right, left)| left.map(|left| PortPair { left, right }))
                .collect::<Vec<_>>()
        };
        let score = |pairs: &[PortPair]| {
            (
                pairs.len(),
                pairs.iter().filter(|&&pair| pair_cost(pair) == 0).count(),
            )
        };
        let unmatched_pairs = best_matching(None);
        // A competing plan must tie on both contraction count and open-pair preference.
        if unmatched_pairs
            .iter()
            .any(|&pair| score(&best_matching(Some(pair))) == score(&unmatched_pairs))
        {
            return Err(TensorCompositionError::Ambiguous(remaining));
        }
        pairs.extend(unmatched_pairs);
        if pairs.is_empty() {
            Ok(ProductPlan::Outer)
        } else {
            Ok(ProductPlan::Contract(pairs))
        }
    }

    /// Select scalar, outer, contraction, or matrix composition from the stored ports.
    /// Ambiguous pairings return an error. The fallible `*` operator delegates here:
    /// `let product = (&left * &right)?;`.
    /// See the [example](SymbolicTensor).
    ///
    /// # Errors
    ///
    /// Rejects invalid source metadata, ambiguous compatible pairings, unsupported
    /// channel composition, or excessive explicit-index multiplicity.
    pub fn multiply(
        &self,
        right: &SymbolicTensor<PartialStructure>,
    ) -> Result<SymbolicTensor<PartialStructure>, TensorCompositionError> {
        for source in [self, right] {
            if source.structure.open_positions().is_empty() {
                source
                    .explicit_index_occurrences()
                    .map_err(TensorCompositionError::from)?;
            } else {
                ExplicitIndexOccurrences::product(&[source])?;
            }
        }
        let plan = self.product_plan(right)?;
        let compact_dot = matches!(&plan, ProductPlan::Contract(_))
            && self.rank() == 1
            && right.rank() == 1
            && [self, right].into_iter().all(|operand| {
                matches!(operand.expression.as_view(), AtomView::Fun(vector)
                    if vector.get_symbol().has_tag(&SPENSO_TAG.rank1))
            });
        if !matches!(&plan, ProductPlan::Compose(..))
            && self.borrowed_parts().product_is_normalized()
            && right.borrowed_parts().product_is_normalized()
        {
            // The plan already checked pairing and ambiguity. Retain its
            // logical order, including nominal ports carried only by zero.
            let interface = match &plan {
                ProductPlan::Scalar if self.is_scalar() => right.structure.clone(),
                ProductPlan::Scalar => self.structure.clone(),
                ProductPlan::Outer => concatenate_interfaces(
                    self.structure.logical_slots(),
                    right.structure.logical_slots(),
                ),
                ProductPlan::Contract(pairs) => concatenate_interfaces(
                    without_positions(
                        &self.structure,
                        &pairs.iter().map(|pair| pair.left).collect::<Vec<_>>(),
                    )
                    .logical_slots(),
                    without_positions(
                        &right.structure,
                        &pairs.iter().map(|pair| pair.right).collect::<Vec<_>>(),
                    )
                    .logical_slots(),
                ),
                ProductPlan::Compose(..) => unreachable!(),
            };
            let occurrences = ExplicitIndexOccurrences::product(&[self, right])?;
            let expression = if compact_dot {
                self.compact_dot_expression(right)?
            } else {
                &self.expression * &right.expression
            };
            let value = Self::from_validated_parts(expression, interface);
            let _ = value.proofs.occurrences.set(Arc::new(occurrences));
            let _ = value.proofs.normalized_product.set(true);
            #[cfg(debug_assertions)]
            if compact_dot {
                TensorProofs::assert_normalization_preserves_interface(
                    "multiply",
                    &[&self.expression, &right.expression],
                    &value.expression,
                    Ok(&value.structure),
                    None,
                );
            }
            return Ok(value);
        }
        plan.apply(self, right)
    }

    /// See the [example](Self::compose).
    ///
    /// # Errors
    ///
    /// Rejects invalid, degenerate, or incompatible channels, unsupported nested
    /// chain scopes, or excessive explicit-index multiplicity.
    pub fn trace_ports(
        &self,
        channel: MatrixChannel,
    ) -> Result<SymbolicTensor<PartialStructure>, TensorCompositionError> {
        if self.structure.open_positions().is_empty() {
            self.explicit_index_occurrences()
                .map_err(TensorCompositionError::from)?;
        } else {
            self.expression.validate_chain_like_nesting()?;
        }
        if channel.input == channel.output {
            return Err(TensorCompositionError::DegenerateChannel {
                input: channel.input,
                output: channel.output,
            });
        }
        let slots = self.structure.logical_slots();
        let input = validate_position(&slots, channel.input)?;
        let output = validate_position(&slots, channel.output)?;
        if !representations_match(&input, &output) {
            return Err(TensorCompositionError::IncompatiblePorts {
                left: channel.input,
                right: channel.output,
            });
        }
        let input_rep = representation(&input);
        let output_rep = representation(&output);
        if !(input_rep.rep.is_self_dual() || input_rep.rep.is_base() && output_rep.rep.is_dual()) {
            return Err(TensorCompositionError::InvalidChannelOrientation {
                input: channel.input,
                output: channel.output,
            });
        }
        let shared = shared_index(
            input,
            output,
            PortPair {
                left: channel.input,
                right: channel.output,
            },
        )?;
        let mut occurrences = ExplicitIndexOccurrences::product(&[self])?;
        if let Some(index) = shared {
            occurrences.add(
                input.rep().slot(index),
                usize::from(matches!(input.aind, PartialIndex::Open(_)))
                    + usize::from(matches!(output.aind, PartialIndex::Open(_))),
            );
            occurrences.validate()?;
        }
        let interface = without_positions(&self.structure, &[channel.input, channel.output]);
        if self.expression.as_view().is_zero() {
            let value = SymbolicTensor::from_normalized_parts(Atom::Zero, interface);
            let _ = value.proofs.occurrences.set(Arc::new(occurrences));
            return value.validate_rewrite("trace_ports", &[self]);
        }

        if !root_chain_channel_is_live(self, channel)
            && matches!(self.expression.as_view(), AtomView::Fun(function)
            if function.get_symbol() == SPENSO_TAG.chain
                || function.get_symbol() == SPENSO_TAG.trace)
        {
            // The root shorthand already owns the global `in`/`out` placeholders.
            // Close a spectator pair in place instead of adding a second consumer.
            let index = shared.unwrap_or_else(|| Self::reserved_dummies([self]).fresh_index());
            let atom = self.materialize_interface_ports(&HashMap::from([
                (channel.input, index),
                (channel.output, index),
            ]))?;
            let interface = if self.structure.open_positions().is_empty() {
                interface
            } else {
                interface_in_atom_order(&atom, &interface)?
            };
            let value = SymbolicTensor::new(atom, interface);
            let _ = value.proofs.occurrences.set(Arc::new(occurrences));
            return value.validate_rewrite("trace_ports", &[self]);
        }
        if contains_chain_like(self.expression.as_view())
            && !root_chain_channel_is_live(self, channel)
        {
            return Err(TensorCompositionError::NestedChainLike);
        }

        let (coefficient, factors) = self.chain_parts(channel)?;
        let atom = coefficient * shadowing::trace(input.rep().to_symbolic([]), factors);
        let interface = if self.structure.open_positions().is_empty() {
            interface
        } else {
            interface_in_atom_order(&atom, &interface)?
        };
        let value = SymbolicTensor::from_normalized_parts(atom, interface);
        let _ = value.proofs.occurrences.set(Arc::new(occurrences));
        value.validate_rewrite("trace_ports", &[self])
    }

    /// See the [example](Self::compose).
    ///
    /// # Errors
    ///
    /// Rejects tensors without a unique matrix channel or a channel that fails tracing.
    pub fn trace_unique(&self) -> Result<SymbolicTensor<PartialStructure>, TensorCompositionError> {
        let channel = self
            .matrix_channel()
            .ok_or(TensorCompositionError::NoMatrixChannel)?;
        self.trace_ports(channel)
    }
}

impl Mul<SymbolicTensor<PartialStructure>> for SymbolicTensor<PartialStructure> {
    type Output = Result<SymbolicTensor<PartialStructure>, TensorCompositionError>;

    #[inline]
    fn mul(self, rhs: SymbolicTensor<PartialStructure>) -> Self::Output {
        self.multiply(&rhs)
    }
}

impl Mul<&SymbolicTensor<PartialStructure>> for SymbolicTensor<PartialStructure> {
    type Output = Result<SymbolicTensor<PartialStructure>, TensorCompositionError>;

    #[inline]
    fn mul(self, rhs: &SymbolicTensor<PartialStructure>) -> Self::Output {
        self.multiply(rhs)
    }
}

impl Mul<SymbolicTensor<PartialStructure>> for &SymbolicTensor<PartialStructure> {
    type Output = Result<SymbolicTensor<PartialStructure>, TensorCompositionError>;

    #[inline]
    fn mul(self, rhs: SymbolicTensor<PartialStructure>) -> Self::Output {
        self.multiply(&rhs)
    }
}

impl Mul<&SymbolicTensor<PartialStructure>> for &SymbolicTensor<PartialStructure> {
    type Output = Result<SymbolicTensor<PartialStructure>, TensorCompositionError>;

    #[inline]
    fn mul(self, rhs: &SymbolicTensor<PartialStructure>) -> Self::Output {
        self.multiply(rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::super::inference::tests::{
        INFERENCE_CALLS, OCCURRENCE_FUNCTION_VISITS, SCOPE_VALIDATIONS,
    };
    use super::*;
    use crate::color::CS;
    use crate::dirac::AGS;
    use crate::representations::{Bispinor, ColorAdjoint, ColorFundamental};
    use crate::test_support;
    use spenso::algebra::ScalarMul;
    use spenso::network::graph::{NetworkNode, NetworkOp};
    use spenso::network::parsing::{NetworkParse, ParseSettings};
    use spenso::structure::dimension::Dimension;
    use spenso::structure::partial::PartialIndex;
    use spenso::structure::representation::{ExtendibleReps, Minkowski, RepName};
    use spenso::structure::slot::{DummyAind, ParseableAind};
    use spenso::{tensor_symbol, vector_symbol};
    use std::convert::Infallible;
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    use std::sync::{Arc, Mutex};
    use symbolica::atom::Symbol;
    use symbolica::{function, symbol};

    fn rep() -> Representation<LibraryRep> {
        ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4))
    }

    fn tensor(name: &str, reps: &[Representation<LibraryRep>]) -> SymbolicTensor<PartialStructure> {
        let symbol: Symbol = SPENSO_TAG.tensor_symbol(name);
        let ports = reps
            .iter()
            .enumerate()
            .map(|(position, rep)| rep.slot(PartialIndex::open(position)))
            .collect::<Vec<_>>();
        partial_tensor(symbol, &ports, &ports)
    }

    fn partial_tensor(
        symbol: Symbol,
        atom_ports: &[PartialSlot],
        logical_ports: &[PartialSlot],
    ) -> SymbolicTensor<PartialStructure> {
        let atom = atom_ports
            .iter()
            .fold(FunctionBuilder::new(symbol), |builder, &slot| {
                builder.add_arg(port_atom(slot))
            })
            .finish();
        SymbolicTensor::new(
            atom,
            PartialStructure::from_logical_slots(logical_ports.iter().copied()),
        )
    }

    fn vector(name: &str, port: PartialSlot) -> SymbolicTensor<PartialStructure> {
        partial_tensor(SPENSO_TAG.rank_one_tensor_symbol(name), &[port], &[port])
    }

    fn factor_slots(atom: &Atom, symbol: Symbol) -> Vec<Slot<LibraryRep, AbstractIndex>> {
        let factors = match atom.as_view() {
            AtomView::Fun(product) if product.get_symbol() == SPENSO_TAG.bracket => {
                product.iter().collect::<Vec<_>>()
            }
            AtomView::Mul(product) => product.iter().collect(),
            _ => panic!("expected a tensor product"),
        };
        factors
            .into_iter()
            .find_map(|factor| {
                let AtomView::Fun(fun) = factor else {
                    return None;
                };
                (fun.get_symbol() == symbol).then(|| {
                    fun.iter()
                        .filter_map(|argument| {
                            Slot::<LibraryRep, AbstractIndex>::try_from(argument).ok()
                        })
                        .collect()
                })
            })
            .unwrap()
    }

    #[test]
    fn bracket_products_flatten_both_sides_without_crossing_other_heads() {
        let a = tensor("flat_a", &[rep()]);
        let b = tensor("flat_b", &[rep()]);
        let c = tensor("flat_c", &[rep()]);
        let d = tensor("flat_d", &[rep()]);
        let left = a.outer(&b).unwrap();
        let right = c.outer(&d).unwrap();
        let result = left.outer(&right).unwrap();
        let expected = &a.expression * &b.expression * &c.expression * &d.expression;
        assert_eq!(result.expression, expected);
        assert_eq!(result.rank(), 4);
        assert_eq!(
            left.outer(&c).unwrap().outer(&d).unwrap().expression,
            expected
        );

        // Even imported, deeply nested brackets flatten, but a sum or a binder
        // remains a single factor (no distribution or placeholder capture).
        let sum = &a.expression + &b.expression;
        let chain = SPENSO_TAG.chain(
            port_atom(a.structure.logical_slots()[0]),
            port_atom(b.structure.logical_slots()[0]),
            [&a.expression, &b.expression],
        );
        let trace = shadowing::trace(rep().to_symbolic([]), [&c.expression, &d.expression]);
        let nested = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(
                FunctionBuilder::new(SPENSO_TAG.bracket)
                    .add_arg(&sum)
                    .add_arg(&chain)
                    .finish(),
            )
            .add_arg(&trace)
            .finish();
        let flattened = SymbolicTensor::bracket_product(nested.as_view(), d.expression.as_view());
        let expected = sum * chain * trace * &d.expression;
        assert_eq!(flattened, expected);
        assert_eq!(
            SymbolicTensor::normalize_products(flattened.as_view()).into_owned(),
            flattened
        );
    }

    #[test]
    fn bracket_removal_preserves_equal_sum_grouping_and_scalar_spectators() {
        let x = Atom::var(symbolica::symbol!("sum_grouping::x"));
        let y = Atom::var(symbolica::symbol!("sum_grouping::y"));
        let s = Atom::var(symbolica::symbol!("sum_grouping::s"));
        let t = Atom::var(symbolica::symbol!("sum_grouping::t"));
        let sum = &x + &y;
        let bracket = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(&sum)
            .finish();
        let spectator = (s + t).pow(8);
        let input = &spectator * (&sum + bracket);
        let expected = Atom::num(2) * &spectator * sum;
        let output = SymbolicTensor::new(input, PartialStructure::from_logical_slots([]));
        assert_eq!(output.expression, expected);
        let retained = SymbolicTensor::from_normalized_parts(
            output.expression.clone(),
            output.structure.clone(),
        );
        assert_eq!(retained, output);
        assert_eq!(
            SymbolicTensor::new(output.expression.clone(), output.structure).expression,
            expected
        );
    }

    #[test]
    fn product_normalization_preserves_opaque_metadata_and_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };

        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = symbolica::symbol!(
            "product_normalization_metadata",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let x = Atom::var(symbolica::symbol!("product_normalization::x"));
        let y = Atom::var(symbolica::symbol!("product_normalization::y"));
        let metadata = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(&x)
            .add_arg(&y)
            .finish();
        let container = FunctionBuilder::new(head).add_arg(metadata).finish();
        let input = (&x + &y).pow(5) * &container + x * y * container;
        calls.store(0, Ordering::Relaxed);

        let normalized = SymbolicTensor::normalize_products(input.as_view());
        assert!(matches!(normalized, AtomOrView::View(_)));
        assert_eq!(normalized.as_view(), input.as_view());
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn product_normalization_flattens_brackets_under_powers_without_merging_open_ports() {
        let vector = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("normalization_open_vector"))
            .add_arg(rep().to_symbolic([]))
            .finish();
        let single = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(&vector)
            .finish();
        let nested = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(single)
            .add_arg(&vector)
            .finish();
        let flat = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(&vector)
            .add_arg(&vector)
            .finish();
        let scalar = Atom::var(symbolica::symbol!("product_normalization::spectator"));
        let input = &scalar * (nested + &scalar).pow(3);
        let expected = &scalar * (flat + &scalar).pow(3);
        let normalized = SymbolicTensor::normalize_products(input.as_view()).into_owned();
        assert_eq!(normalized, expected);
        assert_eq!(
            SymbolicTensor::normalize_products(normalized.as_view()).into_owned(),
            expected
        );
    }

    #[test]
    fn product_normalization_keeps_contracted_power_scope() {
        crate::test_support::test_initialize();
        let slot = rep().slot::<AbstractIndex, _>(73123).to_atom();
        let left = FunctionBuilder::new(spenso::vector_symbol!("inverse_scope_left"))
            .add_arg(&slot)
            .finish();
        let right = FunctionBuilder::new(spenso::vector_symbol!("inverse_scope_right"))
            .add_arg(&slot)
            .finish();
        let product = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(left)
            .add_arg(right)
            .finish();
        for exponent in [-1, -2, -3, 2, 3, 7] {
            let expression = product.pow(exponent);
            let value = SymbolicTensor::from_produced_parts(
                expression.clone(),
                PartialStructure::from_logical_slots([]),
            )
            .unwrap();
            assert_eq!(value.expression, expression);
            assert!(matches!(
                SymbolicTensor::normalize_products(value.expression.as_view()),
                AtomOrView::View(_)
            ));
        }
        // Ordinary scalar brackets retain their existing normalization.
        let x = Atom::var(symbolica::symbol!("inverse_scope_x"));
        let y = Atom::var(symbolica::symbol!("inverse_scope_y"));
        let expression = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(&x)
            .add_arg(&y)
            .finish()
            .pow(-1);
        assert_eq!(
            SymbolicTensor::normalize_products(expression.as_view()).into_owned(),
            (&x * &y).pow(-1)
        );
    }

    #[test]
    fn shared_products_reject_hidden_multiplicity_and_explicit_outer_contractions() {
        let port = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(439)));
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("shared_product_left"),
            &[port],
            &[port],
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("shared_product_right"),
            &[port],
            &[port],
        );
        assert!(matches!(
            left.outer(&right),
            Err(TensorCompositionError::ExplicitOuterContraction(_))
        ));
        let scalar = SymbolicTensor::new(
            &left.expression * &right.expression,
            PartialStructure::from_logical_slots([]),
        );
        let extra = partial_tensor(
            SPENSO_TAG.tensor_symbol("shared_product_extra"),
            &[port],
            &[port],
        );
        for result in [&scalar * &extra, scalar.outer(&extra)] {
            assert!(matches!(
                result,
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            ));
        }
    }

    #[test]
    fn bulk_relabeling_treats_admitted_index_payloads_as_opaque() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let label = symbolica::symbol!(
            "bulk_relabel_index_payload",
            tag = &SPENSO_TAG.index,
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let port = rep().slot(PartialIndex::Explicit(AbstractIndex::Named(
            label.into(),
            7,
            1,
        )));
        let target = rep().slot(PartialIndex::Explicit(AbstractIndex::Named(
            label.into(),
            7,
            2,
        )));
        let head = SPENSO_TAG.tensor_symbol("bulk_relabel_opaque_index_tensor");
        let term = |n, slot| {
            FunctionBuilder::new(head)
                .add_arg(n)
                .add_arg(port_atom(slot))
                .finish()
        };
        let value = SymbolicTensor::infer(Atom::add_many((0..3).map(|n| term(n, port)))).unwrap();
        let expected = Atom::add_many((0..3).map(|n| term(n, target)));
        let replacements = HashMap::from([(0, port_atom(target))]);
        calls.store(0, Ordering::Relaxed);
        assert!(PortRewriteState::exact_arithmetic(
            value.expression.as_view()
        ));
        assert_eq!(
            rewrite_interface_ports(&value, &replacements, None).unwrap(),
            expected
        );
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn exact_bulk_relabeling_matches_ordered_folds_and_checks_every_branch() {
        let port = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(401)));
        let spectator = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(403)));
        let target = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(409)));
        let f = SPENSO_TAG.tensor_symbol("bulk_relabel_f");
        let h = SPENSO_TAG.tensor_symbol("bulk_relabel_h");
        let term = |n| {
            FunctionBuilder::new(f)
                .add_arg(n)
                .add_arg(port_atom(port))
                .finish()
                * FunctionBuilder::new(h)
                    .add_arg(n)
                    .add_arg(port_atom(spectator))
                    .finish()
        };
        let atom = Atom::add_many((1..33).map(term));
        let slots = [port, spectator];
        let replacements = HashMap::from([(0, port_atom(target))]);
        assert!(PortRewriteState::exact_arithmetic(atom.as_view()));
        let mut old = PortRewriteState::new(2);
        let expected = rewrite_ports(atom.as_view(), &slots, &replacements, &mut old).into_owned();
        let value = SymbolicTensor::new(atom.clone(), PartialStructure::from_logical_slots(slots));
        assert_eq!(
            rewrite_interface_ports(&value, &replacements, None).unwrap(),
            expected
        );
        assert_ne!(expected, atom);
        let missing = SymbolicTensor::new(atom + Atom::num(1), value.structure);
        assert!(matches!(
            rewrite_interface_ports(&missing, &replacements, None),
            Err(TensorCompositionError::MissingInterfacePort { position: 0 })
        ));

        let rounded = Atom::num(symbolica::domains::float::Float::parse("0.1", Some(11)).unwrap());
        for atom in [
            rounded.clone() * term(1),
            FunctionBuilder::new(f)
                .add_arg(rounded)
                .add_arg(port_atom(port))
                .finish(),
        ] {
            assert!(!PortRewriteState::exact_arithmetic(atom.as_view()));
            let mut old = PortRewriteState::new(2);
            let expected =
                rewrite_ports(atom.as_view(), &slots, &replacements, &mut old).into_owned();
            let value = SymbolicTensor::new(atom, PartialStructure::from_logical_slots(slots));
            assert_eq!(
                rewrite_interface_ports(&value, &replacements, None).unwrap(),
                expected
            );
        }
    }

    #[test]
    fn relabeling_preserves_callback_interface_and_ranked_zero() {
        let zero = Arc::new(AtomicBool::new(false));
        let emit_zero = Arc::clone(&zero);
        let target = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(419)));
        let target_atom = port_atom(target);
        let head = tensor_symbol!(
            "relabel_result_rank",
            norm = move |function, output| {
                if let AtomView::Fun(function) = function
                    && function.iter().any(|arg| arg == target_atom.as_view())
                {
                    **output = if emit_zero.load(Ordering::Relaxed) {
                        Atom::Zero
                    } else {
                        Atom::num(7)
                            * FunctionBuilder::new(tensor_symbol!("relabel_compliant_result"))
                                .add_args(function.iter())
                                .finish()
                    };
                }
            }
        );
        let port = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(421)));
        let value = partial_tensor(head, &[port], &[port]);
        assert!(!PortRewriteState::exact_arithmetic(
            value.expression.as_view()
        ));
        let replacements = HashMap::from([(0, AbstractIndex::Normal(419))]);
        let rewritten = value.reindex_interface_ports(&replacements).unwrap();
        assert_eq!(
            SymbolicTensor::infer(rewritten.expression.clone())
                .unwrap()
                .structure,
            rewritten.structure
        );
        assert_eq!(rewritten.structure.logical_slots(), [target]);
        assert_eq!(
            rewritten
                .explicit_index_occurrences()
                .unwrap()
                .count(rep().slot(AbstractIndex::Normal(419))),
            1
        );
        let renamed = rewritten
            .reindex_interface_ports(&HashMap::from([(0, AbstractIndex::Normal(425))]))
            .unwrap();
        assert_eq!(
            SymbolicTensor::infer(renamed.expression.clone())
                .unwrap()
                .structure,
            renamed.structure
        );
        assert_eq!(
            renamed.structure.logical_slots(),
            [rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(425)))]
        );
        assert_eq!(
            renamed
                .explicit_index_occurrences()
                .unwrap()
                .count(rep().slot(AbstractIndex::Normal(425))),
            1
        );
        zero.store(true, Ordering::Relaxed);
        let rewritten = value.reindex_interface_ports(&replacements).unwrap();
        assert!(rewritten.expression.as_view().is_zero());
        assert_eq!(rewritten.structure.logical_slots(), [target]);

        let unchanged = tensor_symbol!("relabel_open_pair", norm = |_, _| {});
        let ports = [
            rep().slot(PartialIndex::open(0)),
            rep().slot(PartialIndex::open(1)),
        ];
        let pair = partial_tensor(unchanged, &ports, &ports);
        let paired = pair
            .reindex_interface_ports(&HashMap::from([
                (0, AbstractIndex::Normal(423)),
                (1, AbstractIndex::Normal(423)),
            ]))
            .unwrap();
        assert!(
            SymbolicTensor::from_produced_parts(paired.expression, paired.structure)
                .unwrap()
                .is_scalar()
        );
    }

    #[test]
    fn composition_validation_does_not_replay_materialization_callbacks() {
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let head = tensor_symbol!(
            "composition_observation_callback",
            norm = move |value, _| {
                if let AtomView::Fun(function) = value {
                    observed
                        .lock()
                        .unwrap()
                        .push(function.iter().next().unwrap().to_owned());
                }
            }
        );
        let explicit = PartialIndex::Explicit(AbstractIndex::Normal(431));
        let right_port = rep().slot(explicit);
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("composition_observation_plain"),
            &[right_port],
            &[right_port],
        );
        for left_index in [explicit, PartialIndex::open(0)] {
            let left_port = rep().slot(left_index);
            let left = partial_tensor(head, &[left_port], &[left_port]);
            calls.lock().unwrap().clear();
            let result = if left_index == explicit {
                left.multiply(&right)
            } else {
                left.contract_ports(&right, &[PortPair { left: 0, right: 0 }])
            }
            .unwrap();
            assert!(result.is_scalar());
            let expected = if left_index == explicit {
                Vec::new()
            } else {
                vec![port_atom(right_port)]
            };
            assert_eq!(
                *calls.lock().unwrap(),
                expected,
                "only the actual AUTO port substitution rebuilds the tensor's callback"
            );
        }
    }

    #[test]
    fn scalar_composition_observes_compact_callback_leaves_without_materializing() {
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = tensor_symbol!(
            "scalar_composition_observation_callback",
            norm = move |value, output| {
                observed.fetch_add(1, Ordering::Relaxed);
                if let AtomView::Fun(function) = value
                    && function.iter().any(|argument| {
                        Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok()
                    })
                {
                    **output = Atom::Zero;
                }
            }
        );
        let port = rep().slot(PartialIndex::open(0));
        let tensor = partial_tensor(head, &[port], &[port]);
        let scalar = SymbolicTensor::new(
            Atom::var(symbol!("composition_observation_x")) + Atom::num(1),
            PartialStructure::from_logical_slots([]),
        );
        calls.store(0, Ordering::Relaxed);
        let result = (&scalar * &tensor).unwrap();
        assert_eq!(result.structure, tensor.structure);
        assert!(!result.expression.as_view().is_zero());
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn materialization_and_permutation_keep_occurrence_local_ports() {
        let value = tensor("shared_permutation", &[rep(), rep()]);
        let materialized = value.materialized().unwrap();
        let slots = materialized.structure.logical_slots();
        assert_ne!(slots[0].aind, slots[1].aind);
        assert!(
            slots
                .iter()
                .all(|slot| matches!(slot.aind, PartialIndex::Explicit(_)))
        );
        let reversed = materialized.permuted(&[1, 0]).unwrap();
        assert_eq!(reversed.structure.logical_slots(), [slots[1], slots[0]]);
        assert_eq!(reversed.expression, materialized.expression);
        assert!(matches!(
            value.permuted(&[0, 0]),
            Err(TensorCompositionError::InvalidAxisPermutation { .. })
        ));
        let zero = SymbolicTensor::new(Atom::Zero, value.structure);
        let zero = zero.materialized().unwrap().permuted(&[1, 0]).unwrap();
        assert!(zero.expression.as_view().is_zero());
        assert_eq!(zero.rank(), 2);
    }

    #[test]
    fn identity_port_rewrites_still_validate_every_additive_branch() {
        let port = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(131)));
        let other = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(137)));
        let a = partial_tensor(
            SPENSO_TAG.tensor_symbol("identity_rewrite_a"),
            &[port],
            &[port],
        );
        let b = partial_tensor(
            SPENSO_TAG.tensor_symbol("identity_rewrite_b"),
            &[port],
            &[port],
        );
        let replacements = HashMap::from([(0, port_atom(port))]);
        let atom = &a.expression + &b.expression;
        let mut state = PortRewriteState::new(1);
        let rewritten = rewrite_ports(atom.as_view(), &[port], &replacements, &mut state);
        assert!(matches!(rewritten, AtomOrView::View(_)));
        assert_eq!(rewritten.as_view(), atom.as_view());
        assert_eq!(state.claimed, [true]);
        assert_eq!(state.applied, [true]);

        let missing = partial_tensor(
            SPENSO_TAG.tensor_symbol("identity_rewrite_missing"),
            &[other],
            &[other],
        );
        let invalid = SymbolicTensor::new(&a.expression + missing.expression, a.structure);
        assert!(matches!(
            rewrite_interface_ports(&invalid, &replacements, None),
            Err(TensorCompositionError::MissingInterfacePort { position: 0 })
        ));
    }

    #[test]
    fn port_rewrites_preserve_unchanged_custom_normalizer_calls() {
        let calls = Arc::new(AtomicUsize::new(0));
        let zero = Arc::new(AtomicBool::new(false));
        let observed = Arc::clone(&calls);
        let enabled = Arc::clone(&zero);
        let head = tensor_symbol!(
            "port_rewrite_callback",
            norm = move |_, output| {
                observed.fetch_add(1, Ordering::Relaxed);
                if enabled.load(Ordering::Relaxed) {
                    **output = Atom::Zero;
                }
            }
        );
        let port = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(139)));
        let atom = FunctionBuilder::new(head).add_arg(port_atom(port)).finish();
        for replacements in [HashMap::new(), HashMap::from([(0, port_atom(port))])] {
            calls.store(0, Ordering::Relaxed);
            let mut state = PortRewriteState::new(1);
            let rewritten = rewrite_ports(atom.as_view(), &[port], &replacements, &mut state);
            assert!(matches!(rewritten, AtomOrView::Atom(_)));
            assert_eq!(rewritten.as_view(), atom.as_view());
            assert_eq!(calls.load(Ordering::Relaxed), 1);
            assert_eq!(state.claimed, [true]);
            assert_eq!(state.applied, [replacements.contains_key(&0)]);
        }

        let admitted = SymbolicTensor::from_declared_parts(
            atom.clone(),
            PartialStructure::from_logical_slots([port]),
        )
        .unwrap();
        calls.store(0, Ordering::Relaxed);
        assert_eq!(
            admitted
                .materialize_interface_ports(&HashMap::from([(0, AbstractIndex::Normal(139))]))
                .unwrap(),
            atom
        );
        // Proven identity materialization reuses the admitted payload without
        // rebuilding callback-bearing heads. Raw rewrite scheduling above is
        // still exercised independently.
        assert_eq!(calls.load(Ordering::Relaxed), 0);

        let x = Atom::var(symbol!("port_rewrite_callback_x"));
        let input = &atom * &x + atom.pow(2);
        calls.store(0, Ordering::Relaxed);
        zero.store(true, Ordering::Relaxed);
        let rewritten = rewrite_ports(
            input.as_view(),
            &[],
            &HashMap::new(),
            &mut PortRewriteState::new(0),
        );
        assert!(rewritten.as_view().is_zero());
        assert_eq!(calls.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn port_rewrites_preserve_intrinsic_metric_and_vector_normalization() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let [a, b] = [149, 151]
            .map(|index| rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(index))));
        let metric = FunctionBuilder::new(ETS.metric)
            .add_arg(port_atom(a))
            .add_arg(port_atom(b))
            .finish();
        let x = Atom::var(symbolica::symbol!("port_rewrite_metric_x"));
        let input = (&x + Atom::one()) * &metric;
        let mut state = PortRewriteState::new(2);
        let identity = HashMap::from([(0, port_atom(a)), (1, port_atom(b))]);
        assert!(matches!(
            rewrite_ports(input.as_view(), &[a, b], &identity, &mut state,),
            AtomOrView::View(_)
        ));
        assert_eq!(state.applied, [true, true]);
        let trace = rewrite_ports(
            metric.as_view(),
            &[a, b],
            &HashMap::from([(1, port_atom(a))]),
            &mut PortRewriteState::new(2),
        )
        .into_owned();
        assert_eq!(trace, Atom::num(4));

        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let vector = spenso::vector_symbol!(
            "port_rewrite_vector_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let compact = FunctionBuilder::new(vector)
            .add_arg(mink.to_symbolic([]))
            .finish();
        let dot = FunctionBuilder::new(ETS.metric)
            .add_arg(&compact)
            .add_arg(&compact)
            .finish();
        calls.store(0, Ordering::Relaxed);
        assert!(matches!(
            rewrite_ports(
                dot.as_view(),
                &[],
                &HashMap::new(),
                &mut PortRewriteState::new(0),
            ),
            AtomOrView::View(_)
        ));
        assert_eq!(calls.load(Ordering::Relaxed), 0);

        let replacement = mink.slot(PartialIndex::Explicit(AbstractIndex::Normal(149)));
        let expected = FunctionBuilder::new(vector)
            .add_arg(port_atom(replacement))
            .finish();
        // Reversing source arguments must preserve the mixed metric behavior;
        // symmetry may canonicalize both orders identically. The initially
        // incompatible representations leave the metric intact.
        for arguments in [
            [port_atom(a), compact.clone()],
            [compact.clone(), port_atom(a)],
        ] {
            let mixed = FunctionBuilder::new(ETS.metric)
                .add_args(arguments)
                .finish();
            assert!(
                matches!(mixed.as_view(), AtomView::Fun(fun) if fun.get_symbol() == ETS.metric)
            );
            let unchanged = rewrite_ports(
                mixed.as_view(),
                &[a],
                &HashMap::new(),
                &mut PortRewriteState::new(1),
            );
            assert!(matches!(unchanged, AtomOrView::Atom(_)));
            calls.store(0, Ordering::Relaxed);
            let changed = rewrite_ports(
                mixed.as_view(),
                &[a],
                &HashMap::from([(0, port_atom(replacement))]),
                &mut PortRewriteState::new(1),
            )
            .into_owned();
            assert_eq!(changed, expected);
            assert_eq!(calls.load(Ordering::Relaxed), 1);
        }
    }

    #[test]
    fn port_rewrites_still_normalize_raw_arithmetic_and_functions() {
        let mut sum = Atom::new();
        sum.to_add();
        let mut product = Atom::new();
        product.to_mul();
        let x = Atom::var(symbolica::symbol!("port_rewrite_raw_x"));
        let mut power = Atom::new();
        power.to_pow(x.as_view(), Atom::one().as_view());
        let head = symbolica::symbol!(
            "port_rewrite_raw_function",
            norm = |_, output| {
                **output = Atom::num(7);
            }
        );
        let mut function = Atom::new();
        function.to_fun(head);
        for (input, expected) in [
            (sum, Atom::Zero),
            (product, Atom::one()),
            (power, x),
            (function, Atom::num(7)),
        ] {
            assert!(input.as_view().needs_normalization());
            let rewritten = rewrite_ports(
                input.as_view(),
                &[],
                &HashMap::new(),
                &mut PortRewriteState::new(0),
            )
            .into_owned();
            assert_eq!(rewritten, expected);
            assert!(!rewritten.as_view().needs_normalization());
        }
    }

    #[test]
    fn repeated_explicit_contractions_build_one_nary_product() {
        let slots = (0..6)
            .map(|i| rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(i))))
            .collect::<Vec<_>>();
        let factors = (0..5)
            .map(|i| {
                let ports = &slots[i..i + 2];
                partial_tensor(
                    SPENSO_TAG.tensor_symbol(&format!("nary_factor_{i}")),
                    ports,
                    ports,
                )
            })
            .collect::<Vec<_>>();
        let result = factors[1..].iter().fold(factors[0].clone(), |left, right| {
            left.multiply(right).unwrap()
        });
        let expected = factors
            .iter()
            .fold(Atom::num(1), |product, factor| product * &factor.expression);
        assert_eq!(result.expression, expected);
        assert_eq!(result.structure.logical_slots(), vec![slots[0], slots[5]]);
        let network = result
            .expression
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(
            network
                .graph
                .graph
                .iter_nodes()
                .filter(|(_, _, node)| matches!(node, NetworkNode::Op(NetworkOp::Product)))
                .count(),
            1
        );
        assert_eq!(network.graph.dangling_indices().len(), 2);
    }

    #[test]
    fn admitted_products_reuse_scoped_counts_and_keep_logical_order() {
        test_support::test_initialize();
        let slots = (99100..99106)
            .map(|index| rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(index))))
            .collect::<Vec<_>>();
        let coupling = Atom::var(symbol!("admitted_product_coupling"));
        let factors = (0..5)
            .map(|i| {
                let ports = [slots[i + 1], slots[i]];
                let expression = &coupling
                    * FunctionBuilder::new(
                        SPENSO_TAG.tensor_symbol(&format!("admitted_product_factor_{i}")),
                    )
                    .add_args(ports.map(port_atom))
                    .finish();
                SymbolicTensor::from_declared_parts(
                    expression,
                    PartialStructure::from_logical_slots(ports),
                )
                .unwrap()
            })
            .collect::<Vec<_>>();
        SCOPE_VALIDATIONS.with(|count| count.set(0));
        OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
        assert_eq!(
            factors[0]
                .materialize_interface_ports(&HashMap::from([(0, AbstractIndex::Normal(99101))]))
                .unwrap(),
            factors[0].expression
        );
        let result = factors[1..].iter().fold(factors[0].clone(), |left, right| {
            left.multiply(right).unwrap()
        });
        result.ensure_validated().unwrap();
        assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
        assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
        assert_eq!(
            result.expression,
            factors
                .iter()
                .fold(Atom::one(), |product, factor| product * &factor.expression)
        );
        assert_eq!(result.structure.logical_slots(), [slots[0], slots[5]]);

        // An internal pair remains visible to later products even after its
        // ports have disappeared from the established external interface.
        let extra = SymbolicTensor::infer(
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol("admitted_product_extra"))
                .add_arg(port_atom(slots[2]))
                .finish(),
        )
        .unwrap();
        assert!(matches!(
            result.multiply(&extra),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));
        let produced = SymbolicTensor::from_validated_parts(result.expression, result.structure);
        assert!(matches!(
            produced.multiply(&extra),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));
    }

    #[test]
    fn admitted_product_cancellation_keeps_collision_history_and_local_power_scopes() {
        test_support::test_initialize();
        let scalar = PartialStructure::from_logical_slots([]);
        let leaf = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("admitted_cancelled_power"))
            .add_arg(port_atom(
                rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(99111))),
            ))
            .finish();
        let square = SymbolicTensor::from_declared_parts(leaf.pow(2), scalar.clone()).unwrap();
        let inverse = SymbolicTensor::from_declared_parts(leaf.pow(-2), scalar).unwrap();
        // Construction metadata is authoritative even when ordinary algebra
        // would cancel these factors. Both admitted powers expose two uses of
        // this identity, so their product conservatively rejects the collision.
        assert!(matches!(
            square.multiply(&inverse),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));

        let scalar = SymbolicTensor::from_declared_parts(
            Atom::one() + leaf.pow(2),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        let free = SymbolicTensor::infer(leaf).unwrap();
        let scoped = scalar.multiply(&free).unwrap();
        assert_eq!(scoped.structure, free.structure);
        assert_eq!(
            scoped
                .proofs
                .occurrences
                .get()
                .unwrap()
                .count(rep().slot(AbstractIndex::Normal(99111))),
            1
        );

        // Typed multiplication retains two free oriented ports even if Atom
        // normalization spells their product as a power.
        let fundamental: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let oriented = SymbolicTensor::infer(
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol("admitted_oriented_power"))
                .add_arg(port_atom(
                    fundamental.slot(PartialIndex::Explicit(AbstractIndex::Normal(99113))),
                ))
                .finish(),
        )
        .unwrap();
        let product = oriented.multiply(&oriented).unwrap();
        assert_eq!(product.rank(), 2);
        assert_eq!(
            product.structure.logical_slots(),
            oriented.outer(&oriented).unwrap().structure.logical_slots()
        );
        assert_eq!(product.expression, oriented.expression.pow(2));
    }

    #[test]
    fn repeated_oriented_free_ports_relabel_one_occurrence_and_keep_both_keys() {
        test_support::test_initialize();
        let representation: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let old = AbstractIndex::Normal(99114);
        let new = AbstractIndex::Normal(99115);
        let leaf = |name| {
            SymbolicTensor::infer(
                FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                    .add_arg(representation.slot::<AbstractIndex, _>(old).to_atom())
                    .finish(),
            )
            .unwrap()
        };
        let left = leaf("oriented_free_relabel_v");
        let right = leaf("oriented_free_relabel_w");
        let product = left.multiply(&right).unwrap();
        assert!(product.same_value(&left.outer(&right).unwrap()));
        let relabelled = product
            .reindex_interface_ports(&HashMap::from([(0, new)]))
            .unwrap();
        assert_eq!(
            relabelled.structure.logical_slots(),
            [
                representation.slot(PartialIndex::Explicit(new)),
                representation.slot(PartialIndex::Explicit(old)),
            ]
        );
        let occurrences = relabelled.explicit_index_occurrences().unwrap();
        assert_eq!(occurrences.count(representation.slot(old)), 1);
        assert_eq!(occurrences.count(representation.slot(new)), 1);
        assert_eq!(occurrences.entries().count(), 2);
        assert_eq!(
            occurrences
                .entries()
                .find(|(slot, _)| slot.aind() == old)
                .unwrap()
                .0
                .rep(),
            representation
        );
    }

    #[test]
    fn coalesced_explicit_powers_relabel_independent_complete_copies() {
        test_support::test_initialize();
        let representation: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let [i, k, j, l] = [99161, 99163, 99165, 99167].map(AbstractIndex::Normal);
        let vector = SPENSO_TAG.tensor_symbol("coalesced_copy_vector");
        let atom = |head, indices: &[AbstractIndex]| {
            FunctionBuilder::new(head)
                .add_args(
                    indices
                        .iter()
                        .map(|&index| representation.slot::<AbstractIndex, _>(index).to_atom()),
                )
                .finish()
        };
        let v = SymbolicTensor::infer(atom(vector, &[i])).unwrap();
        let square = v.multiply(&v).unwrap();
        assert!(matches!(square.expression.as_view(), AtomView::Pow(_)));
        for position in [0, 1] {
            let relabelled = square
                .reindex_interface_ports(&HashMap::from([(position, j)]))
                .unwrap();
            assert_eq!(
                relabelled.expression,
                atom(vector, &[j]) * atom(vector, &[i])
            );
            assert_eq!(
                relabelled.structure.logical_slots()[position].aind,
                PartialIndex::Explicit(j)
            );
            assert_eq!(
                relabelled.structure.logical_slots()[1 - position].aind,
                PartialIndex::Explicit(i)
            );
        }
        let matrix = SPENSO_TAG.tensor_symbol("coalesced_copy_multiport");
        let v = SymbolicTensor::infer(atom(matrix, &[i, k])).unwrap();
        let square = v.multiply(&v).unwrap();
        let relabelled = square
            .reindex_interface_ports(&HashMap::from([(0, j), (3, l)]))
            .unwrap();
        assert_eq!(
            relabelled.expression,
            atom(matrix, &[j, k]) * atom(matrix, &[i, l])
        );
        let second = square
            .reindex_interface_ports(&HashMap::from([(2, j)]))
            .unwrap();
        assert_eq!(
            second.expression,
            atom(matrix, &[i, k]) * atom(matrix, &[j, k])
        );
        let permuted = v.permuted(&[1, 0]).unwrap();
        let square = permuted.multiply(&permuted).unwrap();
        let relabelled = square
            .reindex_interface_ports(&HashMap::from([(0, l), (3, j)]))
            .unwrap();
        assert_eq!(
            relabelled.expression,
            atom(matrix, &[i, l]) * atom(matrix, &[j, k])
        );
    }

    #[test]
    fn explicit_relabel_keeps_closed_power_and_additive_indices_local() {
        test_support::test_initialize();
        let i = AbstractIndex::Normal(99171);
        let j = AbstractIndex::Normal(99173);
        let atom = |name, index| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(rep().slot::<AbstractIndex, _>(index).to_atom())
                .finish()
        };
        let a = SymbolicTensor::infer(atom("relabel_closed_a", i)).unwrap();
        let b = SymbolicTensor::infer(atom("relabel_closed_b", i)).unwrap();
        let two = SymbolicTensor::infer(Atom::num(2)).unwrap();
        let power = a.multiply(&b).unwrap().pow(&two).unwrap();
        let sum = SymbolicTensor::infer(Atom::one() + &a.expression * &b.expression).unwrap();
        let v = SymbolicTensor::infer(atom("relabel_closed_v", i)).unwrap();
        for closed in [power, sum] {
            let value = closed.multiply(&v).unwrap();
            let relabelled = value
                .reindex_interface_ports(&HashMap::from([(0, j)]))
                .unwrap();
            assert_eq!(
                relabelled.expression,
                &closed.expression * atom("relabel_closed_v", j)
            );
        }
    }

    #[test]
    fn oriented_branch_bounds_preserve_total_and_hidden_base_after_relabel() {
        test_support::test_initialize();
        let representation: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let i = AbstractIndex::Normal(99181);
        let j = AbstractIndex::Normal(99183);
        let v = SymbolicTensor::infer(
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol("oriented_bounds_v"))
                .add_arg(representation.dual().slot::<AbstractIndex, _>(i).to_atom())
                .finish(),
        )
        .unwrap();
        let base =
            SymbolicTensor::infer(representation.slot::<AbstractIndex, _>(i).to_atom()).unwrap();
        assert!(base.is_scalar());
        let pair = v.multiply(&base).unwrap();
        assert_eq!(pair.rank(), 1);
        let relabelled = pair
            .reindex_interface_ports(&HashMap::from([(0, j)]))
            .unwrap();
        let (hidden, total, base_max, dual_max) = relabelled
            .explicit_index_occurrences()
            .unwrap()
            .oriented_entries()
            .find(|(slot, _, _, _)| slot.aind() == i)
            .unwrap();
        assert_eq!((total, base_max, dual_max), (1, 1, 0));
        assert_eq!(hidden.rep(), representation);
        // Exercise independent branch maxima at their metadata owner. A bare
        // dual slot is tensor syntax, so it is not a scalar-fixture constructor.
        let mut occurrences = ExplicitIndexOccurrences::default();
        occurrences.add(representation.slot(i), 2);
        let mut mixed = ExplicitIndexOccurrences::default();
        mixed.add(representation.slot(i), 1);
        mixed.add(representation.dual().slot(i), 1);
        occurrences.branch_max(&mixed);
        let (_, total, base_max, dual_max) = occurrences.oriented_entries().next().unwrap();
        assert_eq!((total, base_max, dual_max), (2, 2, 1));
        let unchanged = occurrences.mapped::<Infallible>(Ok).unwrap();
        assert_eq!(
            unchanged.oriented_entries().collect::<Vec<_>>(),
            occurrences.oriented_entries().collect::<Vec<_>>()
        );
        let mapped = occurrences
            .mapped::<Infallible>(|slot| {
                Ok(Representation {
                    rep: slot.rep().rep,
                    dim: Dimension::Concrete(if !slot.rep().rep.is_base() { 5 } else { 4 }),
                }
                .slot(slot.aind()))
            })
            .unwrap();
        assert_eq!(mapped.count(representation.rep.new_rep(4).slot(i)), 2);
        assert_eq!(mapped.count(representation.rep.new_rep(5).slot(i)), 1);
    }

    #[test]
    fn admitted_cancellation_and_zero_preserve_hidden_collisions_without_observation() {
        test_support::test_initialize();
        let index = AbstractIndex::Normal(99121);
        let slot = rep().slot(PartialIndex::Explicit(index));
        let leaf = |name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(port_atom(slot))
                .finish()
        };
        let a = SymbolicTensor::infer(leaf("history_cancel_a")).unwrap();
        let b = SymbolicTensor::infer(leaf("history_cancel_b")).unwrap();
        // The inverse closes the already-contracted product in its own scope.
        let minus_one = SymbolicTensor::infer(Atom::num(-1)).unwrap();
        let inverse = a.multiply(&b).unwrap().pow(&minus_one).unwrap();
        let zero = SymbolicTensor::from_declared_parts(
            Atom::Zero,
            PartialStructure::from_logical_slots([slot]),
        )
        .unwrap();
        SCOPE_VALIDATIONS.with(|count| count.set(0));
        OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
        INFERENCE_CALLS.with(|count| count.set(0));

        let pair = a.multiply(&b).unwrap();
        let cancelled = pair.multiply(&inverse).unwrap();
        // A scalar bracket with one argument is an identity. Remove that
        // scope spelling explicitly so ordinary Atom normalization can cancel
        // the ratio, then carry the admitted construction history through it.
        let expression = cancelled.expression.replace_map(|node, _, output| {
            if let AtomView::Fun(function) = node
                && function.get_symbol() == SPENSO_TAG.bracket
                && function.get_nargs() == 1
            {
                **output = function.iter().next().unwrap().to_owned();
            }
        });
        let cancelled = cancelled.with_identity_result(expression, None).unwrap();
        assert_eq!(cancelled.expression, Atom::one());
        assert_eq!(
            cancelled
                .explicit_index_occurrences()
                .unwrap()
                .count(rep().slot(index)),
            2
        );
        assert!(matches!(
            cancelled.multiply(&a),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));

        let zero_pair = zero.multiply(&a).unwrap();
        assert!(zero_pair.is_zero());
        assert!(zero_pair.is_scalar());
        assert_eq!(
            zero_pair
                .explicit_index_occurrences()
                .unwrap()
                .count(rep().slot(index)),
            2
        );
        assert!(matches!(
            zero_pair.multiply(&b),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));
        assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
        assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    }

    #[test]
    fn admitted_compact_dot_retains_explicit_contraction_history() {
        test_support::test_initialize();
        let index = AbstractIndex::Normal(99123);
        let slot = rep().slot(index);
        let vector = |name| {
            SymbolicTensor::infer(
                FunctionBuilder::new(SPENSO_TAG.rank_one_tensor_symbol(name))
                    .add_arg(slot.to_atom())
                    .finish(),
            )
            .unwrap()
        };
        let left = vector("history_dot_p");
        let right = vector("history_dot_q");
        let dot = left.multiply(&right).unwrap();
        assert!(
            matches!(dot.expression.as_view(), AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.dot)
        );
        assert_eq!(dot.explicit_index_occurrences().unwrap().count(slot), 2);
        assert!(matches!(
            dot.multiply(&left),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));
    }

    #[test]
    fn auto_materialization_keeps_cancelled_explicit_history_and_nominal_zero_port() {
        test_support::test_initialize();
        let index = AbstractIndex::Normal(99131);
        let leaf = |name| {
            SymbolicTensor::infer(
                FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                    .add_arg(rep().slot::<AbstractIndex, _>(index).to_atom())
                    .finish(),
            )
            .unwrap()
        };
        let left = leaf("auto_history_left");
        let right = leaf("auto_history_right");
        let zero = SymbolicTensor::infer(Atom::Zero).unwrap();
        let erased = left.multiply(&right).unwrap().multiply(&zero).unwrap();
        let open = SymbolicTensor::infer(
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol("auto_history_open"))
                .add_arg(rep().to_symbolic([]))
                .finish(),
        )
        .unwrap();
        let materialized = erased.outer(&open).unwrap().materialized().unwrap();
        assert!(materialized.is_zero());
        let ports = materialized.structure.logical_slots();
        let [port] = ports.as_slice() else {
            panic!("the nominal zero port survives");
        };
        let PartialIndex::Explicit(materialized_index) = port.aind else {
            panic!("the AUTO port is materialized");
        };
        let occurrences = materialized.explicit_index_occurrences().unwrap();
        assert_eq!(occurrences.count(rep().slot(index)), 2);
        assert_eq!(occurrences.count(rep().slot(materialized_index)), 1);
        assert!(matches!(
            materialized.multiply(&left),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));
    }

    #[test]
    fn occurrence_keys_totally_order_symbolic_dimensions_and_scoped_identities() {
        test_support::test_initialize();
        let index = AbstractIndex::Normal(99125);
        let n = Dimension::Symbolic(symbol!("history_dimension_n").into());
        let m = Dimension::Symbolic(symbol!("history_dimension_m").into());
        let slots = [
            Representation {
                rep: rep().rep,
                dim: n,
            }
            .slot(index),
            rep().slot(index),
            Representation {
                rep: rep().rep,
                dim: m,
            }
            .slot(index),
            Representation {
                rep: rep().rep,
                dim: n,
            }
            .slot(index.scoped(symbol!("history_dimension_scope"))),
        ];
        let mut forward = ExplicitIndexOccurrences::default();
        let mut reverse = ExplicitIndexOccurrences::default();
        for slot in slots {
            forward.add(slot, 1);
        }
        for slot in slots.into_iter().rev() {
            reverse.add(slot, 1);
        }
        assert_eq!(
            forward.entries().collect::<Vec<_>>(),
            reverse.entries().collect::<Vec<_>>()
        );
        assert_eq!(forward.entries().count(), slots.len());
        assert!(
            forward
                .indices
                .windows(2)
                .all(|pair| ExplicitIndexOccurrences::compare(&pair[0], &pair[1]).is_lt())
        );
        for slot in slots {
            assert_eq!(forward.count(slot), 1);
        }
        forward.add(slots[2], 1);
        forward.retain(|slot| *slot != slots[0]);
        assert_eq!(forward.count(slots[0]), 0);
        assert_eq!(forward.count(slots[2]), 2);
        assert_eq!(forward.count(slots[3]), 1);
        let mapped = forward
            .mapped::<Infallible>(|slot| {
                Ok(slot
                    .rep()
                    .slot(slot.aind().scoped(symbol!("history_mapped_scope"))))
            })
            .unwrap();
        assert_eq!(
            mapped.count(
                slots[2]
                    .rep()
                    .slot(slots[2].aind().scoped(symbol!("history_mapped_scope")))
            ),
            2
        );
        mapped.validate().unwrap();
    }

    #[test]
    fn declared_products_still_finish_brackets_and_keep_zero_ports() {
        test_support::test_initialize();
        let x = Atom::var(symbol!("admitted_unfinished_x"));
        let y = Atom::var(symbol!("admitted_unfinished_y"));
        let bracket = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_args([&x, &y])
            .finish();
        let value = SymbolicTensor::from_declared_parts(
            bracket.clone(),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        assert_eq!(value.expression, bracket);
        let one = SymbolicTensor::infer(Atom::one()).unwrap();
        assert_eq!(value.multiply(&one).unwrap().expression, x * y);

        let slot = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(99117)));
        // A declared zero may carry repeated nominal axes. Multiplication by
        // one must not infer a contraction from its absent symbolic payload.
        let zero = SymbolicTensor::from_declared_parts(
            Atom::Zero,
            PartialStructure::from_logical_slots([slot, slot]),
        )
        .unwrap();
        let result = zero.multiply(&one).unwrap();
        assert!(result.is_zero());
        assert_eq!(result.structure.logical_slots(), [slot, slot]);
    }

    #[test]
    fn multiplication_contracts_permuted_explicit_indices_and_preserves_free_ports() {
        let slots = (0..5)
            .map(|index| rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(index))))
            .collect::<Vec<_>>();
        let left_ports = &slots[..4];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("unique_named_left"),
            left_ports,
            left_ports,
        );
        for last in [0, 4] {
            let right_ports = [slots[2], slots[last], slots[3], slots[1]];
            let right = partial_tensor(
                SPENSO_TAG.tensor_symbol("unique_named_right"),
                &right_ports,
                &right_ports,
            );
            let result = left.multiply(&right).unwrap();
            let expected = if last == 0 {
                Vec::new()
            } else {
                vec![slots[0], slots[4]]
            };
            assert_eq!(result.structure.logical_slots(), expected);
            let network = result
                .expression
                .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
                .unwrap();
            assert_eq!(network.graph.dangling_indices().len(), expected.len());
            let zero = SymbolicTensor::new(Atom::Zero, left.structure.clone());
            let zero_result = zero.multiply(&right).unwrap();
            assert!(zero_result.expression.as_view().is_zero());
            assert_eq!(zero_result.structure.logical_slots(), expected);
        }
    }

    #[test]
    fn named_pairs_take_precedence_over_open_ports_and_matrix_channels() {
        let a = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(101)));
        let b = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(102)));
        let c = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(103)));
        let open = rep().slot(PartialIndex::open(0));
        for (left_ports, right_ports) in [
            (vec![a, open], vec![a, b]),
            (vec![a, open, c], vec![a, b, c]),
        ] {
            let left = partial_tensor(
                SPENSO_TAG.tensor_symbol("named_priority_left"),
                &left_ports,
                &left_ports,
            );
            let right = partial_tensor(
                SPENSO_TAG.tensor_symbol("named_priority_right"),
                &right_ports,
                &right_ports,
            );
            let result = left.multiply(&right).unwrap();
            assert!(result.is_scalar());
            assert!(
                result
                    .expression
                    .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
                    .unwrap()
                    .graph
                    .dangling_indices()
                    .is_empty()
            );
        }
    }

    #[test]
    fn overlapping_candidates_with_a_unique_maximum_matching_are_not_ambiguous() {
        let a = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(111)));
        let b = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(112)));
        let c = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(113)));
        let open = rep().slot(PartialIndex::open(0));
        let left_ports = [a, open, c];
        let right_ports = [b, open, c];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("unique_maximum_left"),
            &left_ports,
            &left_ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("unique_maximum_right"),
            &right_ports,
            &right_ports,
        );
        let result = left.multiply(&right).unwrap();
        assert!(result.is_scalar());
        assert!(
            result
                .expression
                .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
                .unwrap()
                .graph
                .dangling_indices()
                .is_empty()
        );
    }

    #[test]
    fn unresolved_pairs_take_precedence_over_unmatched_named_ports() {
        let named = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(114)));
        let open = rep().slot(PartialIndex::open(0));
        let vector = partial_tensor(
            SPENSO_TAG.tensor_symbol("preferred_unresolved_vector"),
            &[open],
            &[open],
        );
        for ports in [[named, open], [open, named]] {
            let matrix = partial_tensor(
                SPENSO_TAG.tensor_symbol("preferred_unresolved_matrix"),
                &ports,
                &ports,
            );
            let open_position = ports.iter().position(|port| *port == open).unwrap();
            for (left, right, expected) in [
                (
                    &vector,
                    &matrix,
                    PortPair {
                        left: 0,
                        right: open_position,
                    },
                ),
                (
                    &matrix,
                    &vector,
                    PortPair {
                        left: open_position,
                        right: 0,
                    },
                ),
            ] {
                assert!(matches!(
                    left.product_plan(right),
                    Ok(ProductPlan::Contract(pairs)) if pairs == [expected]
                ));
                let result = left.multiply(right).unwrap();
                assert_eq!(result.structure.logical_slots(), vec![named]);
                assert_eq!(
                    result
                        .expression
                        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
                        .unwrap()
                        .graph
                        .dangling_indices()
                        .len(),
                    1
                );
            }
        }
    }

    #[test]
    fn unresolved_preference_preserves_ambiguity_between_open_partners() {
        let named = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(115)));
        let open = rep().slot(PartialIndex::open(0));
        let vector = partial_tensor(
            SPENSO_TAG.tensor_symbol("equally_preferred_vector"),
            &[open],
            &[open],
        );
        let ports = [named, open, rep().slot(PartialIndex::open(1))];
        let tensor = partial_tensor(
            SPENSO_TAG.tensor_symbol("equally_preferred_tensor"),
            &ports,
            &ports,
        );
        for (left, right) in [(&vector, &tensor), (&tensor, &vector)] {
            assert!(matches!(
                left.multiply(right),
                Err(TensorCompositionError::Ambiguous(_))
            ));
        }
    }

    #[test]
    fn unresolved_pairing_matches_exhaustive_rank_three_plans() {
        fn enumerate(
            position: usize,
            candidates: &[PortPair],
            pairs: &mut Vec<PortPair>,
            plans: &mut Vec<Vec<PortPair>>,
        ) {
            if position == 3 {
                plans.push(pairs.clone());
                return;
            }
            enumerate(position + 1, candidates, pairs, plans);
            for &pair in candidates.iter().filter(|pair| pair.left == position) {
                if pairs.iter().any(|existing| existing.right == pair.right) {
                    continue;
                }
                pairs.push(pair);
                enumerate(position + 1, candidates, pairs, plans);
                pairs.pop();
            }
        }

        let ports = |mask, offset| {
            (0..3)
                .map(|position| {
                    let index = if mask & (1 << position) != 0 {
                        PartialIndex::open(position)
                    } else {
                        PartialIndex::Explicit(AbstractIndex::Normal(position + offset))
                    };
                    rep().slot(index)
                })
                .collect::<Vec<_>>()
        };
        for left_mask in 0..8 {
            for right_mask in 0..8 {
                for offset in [0, 3] {
                    let left_ports = ports(left_mask, 0);
                    let right_ports = ports(right_mask, offset);
                    let left = partial_tensor(
                        SPENSO_TAG.tensor_symbol("exhaustive_preference_left"),
                        &left_ports,
                        &left_ports,
                    );
                    let right = partial_tensor(
                        SPENSO_TAG.tensor_symbol("exhaustive_preference_right"),
                        &right_ports,
                        &right_ports,
                    );
                    let candidates = compatible_pairs(&left.structure, &right.structure);
                    let named = candidates
                        .iter()
                        .filter(|pair| {
                            matches!(
                                (left_ports[pair.left].aind, right_ports[pair.right].aind),
                                (PartialIndex::Explicit(_), PartialIndex::Explicit(_))
                            )
                        })
                        .collect::<Vec<_>>();
                    let mut plans = Vec::new();
                    enumerate(0, &candidates, &mut Vec::new(), &mut plans);
                    plans.retain(|plan| named.iter().all(|pair| plan.contains(pair)));
                    let score = |pairs: &[PortPair]| {
                        (
                            pairs.len(),
                            pairs
                                .iter()
                                .filter(|pair| {
                                    matches!(
                                        (left_ports[pair.left].aind, right_ports[pair.right].aind),
                                        (PartialIndex::Open(_), PartialIndex::Open(_))
                                    )
                                })
                                .count(),
                        )
                    };
                    let best = plans.iter().map(|plan| score(plan)).max().unwrap();
                    let preferred = plans
                        .iter()
                        .filter(|plan| score(plan) == best)
                        .collect::<Vec<_>>();
                    match left.product_plan(&right) {
                        Err(TensorCompositionError::Ambiguous(_)) => assert!(preferred.len() > 1),
                        Ok(ProductPlan::Outer) => assert_eq!(best, (0, 0)),
                        Ok(ProductPlan::Contract(actual)) => {
                            assert_eq!(preferred.len(), 1);
                            assert_eq!(actual.len(), preferred[0].len());
                            assert!(actual.iter().all(|pair| preferred[0].contains(pair)));
                        }
                        plan => panic!("unexpected plan {plan:?}"),
                    }
                }
            }
        }
    }

    #[test]
    fn distinct_representations_determine_multiple_open_contractions() {
        let reps = [
            rep(),
            Minkowski {}.new_rep(4).cast(),
            ColorAdjoint {}.new_rep(8).cast(),
        ];
        let left = tensor("unique_open_left", &reps);
        let right = tensor("unique_open_right", &[reps[2], reps[0], reps[1]]);
        let result = (left * &right).unwrap();
        assert!(result.is_scalar());
        assert!(
            result
                .expression
                .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
                .unwrap()
                .graph
                .dangling_indices()
                .is_empty()
        );
    }

    #[test]
    fn equally_large_open_matchings_remain_ambiguous_after_named_contraction() {
        let ports = [
            rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(121))),
            rep().slot(PartialIndex::open(0)),
            rep().slot(PartialIndex::open(1)),
        ];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("ambiguous_open_left"),
            &ports,
            &ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("ambiguous_open_right"),
            &ports,
            &ports,
        );
        assert!(
            matches!(&left * right, Err(TensorCompositionError::Ambiguous(pairs)) if pairs.len() == 4)
        );
    }

    #[test]
    fn matrix_multiplication_keeps_open_spectators() {
        let reps = [rep(), rep(), Minkowski {}.new_rep(4).cast()];
        let left = tensor("open_spectator_matrix_left", &reps);
        let right = tensor("open_spectator_matrix_right", &reps);
        let result = (left * right).unwrap();
        assert_eq!(result.rank(), 4);
        assert!(
            matches!(result.expression.as_view(), AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.chain)
        );
    }

    #[test]
    fn automatic_dot_requires_tagged_vector_heads() {
        let port = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(68101)));
        let tagged = vector("automatic_dot_tagged", port);
        let other_tagged = vector("automatic_dot_other_tagged", port);
        let untagged = partial_tensor(
            SPENSO_TAG.tensor_symbol("automatic_dot_untagged"),
            &[port],
            &[port],
        );
        let other_untagged = partial_tensor(
            SPENSO_TAG.tensor_symbol("automatic_dot_other_untagged"),
            &[port],
            &[port],
        );
        for (left, right, uses_dot) in [
            (&tagged, &other_tagged, true),
            (&tagged, &untagged, false),
            (&untagged, &tagged, false),
            (&untagged, &other_untagged, false),
        ] {
            let result = left.multiply(right).unwrap();
            assert!(result.is_scalar());
            assert_eq!(result.expression.contains_symbol(SPENSO_TAG.dot), uses_dot);
            let network = result
                .expression
                .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
                .unwrap();
            assert!(network.state.is_scalar());
            assert!(network.graph.dangling_indices().is_empty());
            let explicit = left.dot(right).unwrap();
            assert!(explicit.expression.contains_symbol(SPENSO_TAG.dot));
            let network = explicit
                .expression
                .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
                .unwrap();
            assert!(network.state.is_scalar());
        }
    }

    #[test]
    fn metric_vector_multiplication_orders_keep_schoonschip_operands() {
        let rep: Representation<LibraryRep> = crate::test_support::test_initialize().mink4.cast();
        let mu = rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(68102)));
        let nu = rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(68103)));
        let p = vector("metric_order_tagged_p", mu);
        let q = vector("metric_order_tagged_q", nu);
        let metric = SymbolicTensor::infer(spenso::g!(port_atom(mu), port_atom(nu))).unwrap();
        let compact_vector = |value: &SymbolicTensor<PartialStructure>| {
            rewrite_interface_ports(value, &HashMap::from([(0, rep.to_symbolic([]))]), None)
                .unwrap()
        };
        let expected = spenso::g!(compact_vector(&p), compact_vector(&q));
        let factors = [&metric, &p, &q];
        let mut common = None;
        for [a, b, c] in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            let product = factors[a]
                .multiply(factors[b])
                .unwrap()
                .multiply(factors[c])
                .unwrap();
            assert!(product.is_scalar());
            assert!(!product.expression.contains_symbol(SPENSO_TAG.dot));
            if let Some(common) = &common {
                assert_eq!(&product, common);
            } else {
                common = Some(product.clone());
            }
            for metrics in [false, true] {
                let settings = crate::tensor::ContractSettings {
                    metrics,
                    ..Default::default()
                };
                let result = product.contract(settings).unwrap();
                assert_eq!(result.expression, expected);
                assert!(result.is_scalar());
                assert_eq!(result.contract(settings).unwrap(), result);
            }
        }
    }

    #[test]
    fn rank_one_product_uses_dot() {
        let left = vector("rank_one_product_left", rep().slot(PartialIndex::open(0)));
        let right = vector("rank_one_product_right", rep().slot(PartialIndex::open(0)));
        let result = left.multiply(&right).unwrap();

        assert!(result.is_scalar());
        assert!(
            matches!(result.expression.as_view(), AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.dot)
        );
    }

    #[test]
    fn explicit_rank_one_contraction_keeps_compact_dot_syntax() {
        let index = AbstractIndex::Normal(7);
        let port = rep().slot(PartialIndex::Explicit(index));
        let left_symbol = SPENSO_TAG.rank_one_tensor_symbol("explicit_dot_left");
        let right_symbol = SPENSO_TAG.rank_one_tensor_symbol("explicit_dot_right");
        let left = partial_tensor(left_symbol, &[port], &[port]);
        let right = partial_tensor(right_symbol, &[port], &[port]);
        let result =
            SymbolicTensor::contract_ports(&left, &right, &[PortPair { left: 0, right: 0 }])
                .unwrap();

        let AtomView::Fun(dot) = result.expression.as_view() else {
            panic!("expected a dot")
        };
        assert_eq!(dot.get_symbol(), SPENSO_TAG.dot);
        for operand in dot.iter() {
            let AtomView::Fun(vector) = operand else {
                panic!("expected a compact vector")
            };
            let argument = vector.iter().next().unwrap();
            assert!(Representation::<LibraryRep>::try_from(argument).is_ok());
            assert!(Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_err());
        }

        let network = result
            .expression
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(network.state.is_scalar());
        assert!(network.graph.dangling_indices().is_empty());
    }

    #[test]
    fn dual_rank_one_dot_materializes_with_opposite_orientations() {
        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let dual = base.dual();
        let left = vector("dual_dot_left", base.slot(PartialIndex::open(0)));
        let right = vector("dual_dot_right", dual.slot(PartialIndex::open(0)));
        let result = left.multiply(&right).unwrap();

        assert!(
            matches!(result.expression.as_view(), AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.dot)
        );
        let network = result
            .expression
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(
            network.state.is_scalar(),
            "state {:?}, dangling {:?}",
            network.state,
            network.graph.dangling_indices()
        );
        assert!(network.graph.dangling_indices().is_empty());
    }

    #[test]
    fn matrix_products_flatten_chains() {
        let first = tensor("first", &[rep(), rep()]);
        let second = tensor("second", &[rep(), rep()]);
        let third = tensor("third", &[rep(), rep()]);
        let chain = first.multiply(&second).unwrap();
        let chain = chain.multiply(&third).unwrap();

        let AtomView::Fun(fun) = chain.expression.as_view() else {
            panic!("expected a chain")
        };
        assert_eq!(fun.get_symbol(), SPENSO_TAG.chain);
        assert_eq!(fun.get_nargs(), 5);
        assert_eq!(chain.rank(), 2);
    }

    #[test]
    fn scalar_chain_wrappers_compose_and_trace_but_broadcasts_remain_scoped() {
        let first = tensor("wrapped_chain_first", &[rep(), rep()]);
        let second = tensor("wrapped_chain_second", &[rep(), rep()]);
        let third = tensor("wrapped_chain_third", &[rep(), rep()]);
        let chain = first.multiply(&second).unwrap();
        let scalar = Atom::var(symbolica::symbol!("wrapped_chain_scalar"));
        let broadcast = FunctionBuilder::new(spenso::broadcast_symbol!("wrapped_chain_broadcast"))
            .add_arg(&chain.expression)
            .finish();

        for (atom, can_extend) in [
            (scalar.clone() * chain.expression.as_ref(), true),
            (broadcast, false),
        ] {
            let wrapped = SymbolicTensor::new(atom, chain.structure.clone());
            assert_eq!(
                SymbolicTensor::matrix_channel(&wrapped),
                Some(MatrixChannel {
                    input: 0,
                    output: 1,
                })
            );

            let product = wrapped.multiply(&third);
            let traced = SymbolicTensor::trace_unique(&wrapped);
            if can_extend {
                let expected = chain.multiply(&third).unwrap();
                assert_eq!(
                    product.unwrap().expression,
                    scalar.clone() * expected.expression
                );
                let expected = SymbolicTensor::trace_unique(&chain).unwrap();
                assert_eq!(
                    traced.unwrap().expression,
                    scalar.clone() * expected.expression
                );
            } else {
                assert!(matches!(
                    product,
                    Err(TensorCompositionError::NestedChainLike)
                ));
                assert!(matches!(
                    traced,
                    Err(TensorCompositionError::NestedChainLike)
                ));
            }
        }
    }

    #[test]
    fn contracted_chain_endpoints_are_not_treated_as_the_public_channel() {
        let channel = rep();
        let endpoint = channel
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(61))
            .to_atom();
        let spectators = [
            ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4)),
            ColorAdjoint {}.new_rep(8).cast(),
            Bispinor {}.new_rep(4).cast(),
        ];
        let factor = spectators.iter().fold(
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol("closed_chain_spectators"))
                .add_arg(Atom::var(SPENSO_TAG.chain_in))
                .add_arg(Atom::var(SPENSO_TAG.chain_out)),
            |builder, representation| builder.add_arg(representation.to_symbolic([])),
        );
        let value = SymbolicTensor::new(
            SPENSO_TAG.chain(&endpoint, &endpoint, [factor.finish()]),
            PartialStructure::from_logical_slots(spectators.into_iter().enumerate().map(
                |(position, representation)| representation.slot(PartialIndex::open(position)),
            )),
        );

        assert_eq!(SymbolicTensor::matrix_channel(&value), None);
    }

    #[test]
    fn matrix_composition_contracts_equal_explicit_spectators() {
        let channel = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(47);
        let left_ports = [
            channel.slot(PartialIndex::open(0)),
            channel.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let right_ports = [
            channel.slot(PartialIndex::open(0)),
            channel.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("explicit_spectator_left"),
            &left_ports,
            &left_ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("explicit_spectator_right"),
            &right_ports,
            &right_ports,
        );

        let result = left.multiply(&right).unwrap();

        assert_eq!(result.rank(), 2);
        assert!(
            matches!(result.expression.as_view(), AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.chain)
        );
        let materialized = SymbolicTensor::materialize_interface_ports(
            &result,
            &HashMap::from([
                (0, AbstractIndex::Normal(61)),
                (1, AbstractIndex::Normal(67)),
            ]),
        )
        .unwrap();
        let network = materialized
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(network.graph.dangling_indices().len(), 2);
    }

    #[test]
    fn explicit_contraction_contracts_equal_explicit_spectators() {
        let channel = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(51);
        let left_ports = [
            channel.slot(PartialIndex::open(0)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let right_ports = [
            channel.slot(PartialIndex::open(0)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("explicit_contract_spectator_left"),
            &left_ports,
            &left_ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("explicit_contract_spectator_right"),
            &right_ports,
            &right_ports,
        );

        let result =
            SymbolicTensor::contract_ports(&left, &right, &[PortPair { left: 0, right: 0 }])
                .unwrap();

        assert!(result.is_scalar());
        let network = result
            .expression
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(network.graph.dangling_indices().is_empty());
    }

    #[test]
    fn matrix_composition_rejects_ambiguous_explicit_spectators() {
        let channel = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(53);
        let left_ports = [
            channel.slot(PartialIndex::open(0)),
            channel.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::Explicit(index)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let right_ports = [
            channel.slot(PartialIndex::open(0)),
            channel.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("ambiguous_spectator_left"),
            &left_ports,
            &left_ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("ambiguous_spectator_right"),
            &right_ports,
            &right_ports,
        );

        assert!(matches!(
            SymbolicTensor::compose(
                &left,
                &right,
                MatrixChannel {
                    input: 0,
                    output: 1,
                },
                MatrixChannel {
                    input: 0,
                    output: 1,
                },
            ),
            Err(TensorCompositionError::Ambiguous(candidates)) if candidates.len() == 2
        ));
    }

    #[test]
    fn multiple_compatible_spectators_are_ambiguous() {
        let left = tensor("ambiguous_spectator_left", &[rep(), rep(), rep()]);
        let right = tensor("ambiguous_spectator_right", &[rep()]);
        let error = left.multiply(&right).unwrap_err();

        assert!(
            matches!(error, TensorCompositionError::Ambiguous(candidates) if candidates.len() == 3)
        );
    }

    #[test]
    fn outer_materialization_uses_logical_operand_order() {
        let q = tensor("z_outer_q", &[rep()]);
        let p = tensor("a_outer_p", &[rep()]);
        let product = q.outer(&p).unwrap();
        let indexed = SymbolicTensor::materialize_interface_ports(
            &product,
            &HashMap::from([(0, AbstractIndex::Normal(17))]),
        )
        .unwrap();
        let q_symbol = SPENSO_TAG.tensor_symbol("z_outer_q");
        let p_symbol = SPENSO_TAG.tensor_symbol("a_outer_p");
        let factors = match indexed.as_view() {
            AtomView::Fun(product) if product.get_symbol() == SPENSO_TAG.bracket => {
                product.iter().collect::<Vec<_>>()
            }
            AtomView::Mul(product) => product.iter().collect(),
            _ => panic!("expected a tensor product"),
        };
        let mut q_is_indexed = false;
        let mut p_is_open = false;
        for factor in factors {
            let AtomView::Fun(fun) = factor else {
                continue;
            };
            let argument = fun.iter().next().unwrap();
            if fun.get_symbol() == q_symbol {
                q_is_indexed = Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok();
            } else if fun.get_symbol() == p_symbol {
                p_is_open = Representation::<LibraryRep>::try_from(argument).is_ok()
                    && Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_err();
            }
        }
        assert!(q_is_indexed && p_is_open);
    }

    #[test]
    fn materialization_matches_logical_ports_in_canonical_atom_order() {
        let euc = rep();
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let euc_port = euc.slot(PartialIndex::open(0));
        let mink_port = mink.slot(PartialIndex::open(1));
        let value = partial_tensor(
            SPENSO_TAG.tensor_symbol("canonical_port_order"),
            &[euc_port, mink_port],
            &[mink_port, euc_port],
        );
        let indexed = SymbolicTensor::materialize_interface_ports(
            &value,
            &HashMap::from([(0, AbstractIndex::Normal(19))]),
        )
        .unwrap();

        let AtomView::Fun(tensor) = indexed.as_view() else {
            panic!("expected a tensor")
        };
        let args = tensor.iter().collect::<Vec<_>>();
        assert_eq!(
            Representation::<LibraryRep>::try_from(args[0]).unwrap(),
            euc
        );
        assert_eq!(
            Slot::<LibraryRep, AbstractIndex>::try_from(args[1]).unwrap(),
            mink.slot(AbstractIndex::Normal(19))
        );
    }

    #[test]
    fn materialization_does_not_rewrite_tensor_scalar_metadata() {
        let representation = rep();
        let port = representation.slot(PartialIndex::open(0));
        let metadata = FunctionBuilder::new(symbolica::symbol!("composition_scalar_metadata"))
            .add_arg(representation.to_symbolic([]))
            .finish();
        let atom = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("tensor_with_scalar_metadata"))
            .add_arg(metadata)
            .add_arg(representation.to_symbolic([]))
            .finish();
        let value = SymbolicTensor::new(atom, PartialStructure::from_logical_slots([port]));

        let indexed = SymbolicTensor::materialize_interface_ports(
            &value,
            &HashMap::from([(0, AbstractIndex::Normal(21))]),
        )
        .unwrap();
        let AtomView::Fun(tensor) = indexed.as_view() else {
            panic!("expected a tensor")
        };
        let arguments = tensor.iter().collect::<Vec<_>>();
        let AtomView::Fun(metadata) = arguments[0] else {
            panic!("expected scalar metadata")
        };

        assert!(Representation::<LibraryRep>::try_from(metadata.iter().next().unwrap()).is_ok());
        assert_eq!(
            Slot::<LibraryRep, AbstractIndex>::try_from(arguments[1]).unwrap(),
            representation.slot(AbstractIndex::Normal(21))
        );
    }

    #[test]
    fn registered_and_imported_projectors_relabel_nested_ports() {
        let a = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(191)));
        let b = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(192)));
        let head = SPENSO_TAG.tensor_symbol("classifier_projector_operand");
        let leaf = partial_tensor(head, &[a], &[a]);
        let renamed = partial_tensor(head, &[b], &[b]);
        for projector in [*shadowing::SYM, *shadowing::ANTISYM, *shadowing::CYCLIC] {
            let expression = projector.call(&leaf.expression);
            let mut bytes = Vec::new();
            expression.export(&mut bytes).unwrap();
            let imported = Atom::import(&mut bytes.as_slice(), None).unwrap();
            assert_eq!(imported, expression);
            for expression in [expression, imported] {
                let value = SymbolicTensor::<PartialStructure>::infer(expression).unwrap();
                let result = value
                    .reindex_interface_ports(&HashMap::from([(0, AbstractIndex::Normal(192))]))
                    .unwrap();
                assert_eq!(result.expression, projector.call(&renamed.expression));
                assert_eq!(result.structure.logical_slots(), vec![b]);
            }
        }
    }

    #[test]
    fn attributed_tensor_leaves_relabel_only_declared_ports() {
        let a = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(193)));
        let b = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(194)));
        let metadata = Atom::var(symbolica::symbol!("classifier_metadata"));
        let heads = [
            SPENSO_TAG.tensor_symbol("classifier_plain_tensor"),
            symbolica::symbol!("classifier_symmetric_tensor"; Symmetric; tags = [&SPENSO_TAG.tensor]),
            symbolica::symbol!("classifier_antisymmetric_tensor"; Antisymmetric; tags = [&SPENSO_TAG.tensor]),
            symbolica::symbol!("classifier_cyclic_tensor"; Cyclesymmetric; tags = [&SPENSO_TAG.tensor]),
        ];
        for head in heads {
            // Matching a projector's attribute does not make this a projector.
            assert!(is_tensor_leaf_head(head));
            let expression = head.call((&metadata, port_atom(a)));
            let value = SymbolicTensor::<PartialStructure>::infer(expression).unwrap();
            let result = value
                .reindex_interface_ports(&HashMap::from([(0, AbstractIndex::Normal(194))]))
                .unwrap();
            assert_eq!(result.expression, head.call((&metadata, port_atom(b))));
            assert_eq!(result.structure.logical_slots(), vec![b]);
        }
    }

    #[test]
    fn materialization_does_not_rewrite_scalar_factor_metadata() {
        let representation = rep();
        let port = representation.slot(PartialIndex::open(0));
        let metadata_symbol = symbolica::symbol!("composition_external_scalar_metadata");
        let tensor_symbol = SPENSO_TAG.tensor_symbol("tensor_beside_scalar_metadata");
        let metadata = FunctionBuilder::new(metadata_symbol)
            .add_arg(representation.to_symbolic([]))
            .finish();
        let tensor = FunctionBuilder::new(tensor_symbol)
            .add_arg(representation.to_symbolic([]))
            .finish();
        let value = SymbolicTensor::new(
            metadata * tensor,
            PartialStructure::from_logical_slots([port]),
        );

        let indexed = SymbolicTensor::materialize_interface_ports(
            &value,
            &HashMap::from([(0, AbstractIndex::Normal(22))]),
        )
        .unwrap();
        let AtomView::Mul(product) = indexed.as_view() else {
            panic!("expected a product")
        };
        let arguments = product
            .iter()
            .filter_map(|factor| {
                let AtomView::Fun(function) = factor else {
                    return None;
                };
                Some((function.get_symbol(), function.iter().next().unwrap()))
            })
            .collect::<HashMap<_, _>>();

        assert!(Representation::<LibraryRep>::try_from(arguments[&metadata_symbol]).is_ok());
        assert_eq!(
            Slot::<LibraryRep, AbstractIndex>::try_from(arguments[&tensor_symbol]).unwrap(),
            representation.slot(AbstractIndex::Normal(22))
        );
    }

    #[test]
    fn structured_zero_preserves_interface_metadata_through_tensor_operations() {
        let representation = rep();
        let vector_port = representation.slot(PartialIndex::open(0));
        let zero_vector = SymbolicTensor::new(
            Atom::Zero,
            PartialStructure::from_logical_slots([vector_port]),
        );
        let vector = tensor("structured_zero_vector_partner", &[representation]);

        assert!(
            SymbolicTensor::materialize_interface_ports(
                &zero_vector,
                &HashMap::from([(0, AbstractIndex::Normal(24))]),
            )
            .unwrap()
            .as_view()
            .is_zero()
        );
        let dotted = SymbolicTensor::contract_ports(
            &zero_vector,
            &vector,
            &[PortPair { left: 0, right: 0 }],
        )
        .unwrap();
        assert!(dotted.is_scalar());
        assert!(dotted.expression.as_view().is_zero());

        let outer_product = zero_vector.outer(&vector).unwrap();
        assert_eq!(outer_product.rank(), 2);
        assert!(outer_product.expression.as_view().is_zero());

        let matrix_ports = [
            representation.slot(PartialIndex::open(0)),
            representation.slot(PartialIndex::open(1)),
        ];
        let zero_matrix = SymbolicTensor::new(
            Atom::Zero,
            PartialStructure::from_logical_slots(matrix_ports),
        );
        let matrix = tensor("structured_zero_matrix_partner", &[representation; 2]);
        let composed = SymbolicTensor::compose(
            &zero_matrix,
            &matrix,
            MatrixChannel {
                input: 0,
                output: 1,
            },
            MatrixChannel {
                input: 0,
                output: 1,
            },
        )
        .unwrap();
        assert_eq!(composed.rank(), 2);
        assert!(composed.expression.as_view().is_zero());

        let traced = SymbolicTensor::trace_ports(
            &zero_matrix,
            MatrixChannel {
                input: 0,
                output: 1,
            },
        )
        .unwrap();
        assert!(traced.is_scalar());
        assert!(traced.expression.as_view().is_zero());
    }

    #[test]
    fn surviving_ports_skip_internal_contraction_dummies() {
        let left_symbol = SPENSO_TAG.tensor_symbol("dummy_skip_left");
        let right_symbol = SPENSO_TAG.tensor_symbol("dummy_skip_right");
        let ports = [
            rep().slot(PartialIndex::open(0)),
            rep().slot(PartialIndex::open(1)),
        ];
        let left = partial_tensor(left_symbol, &ports, &ports);
        let right = partial_tensor(right_symbol, &ports, &ports);
        let contracted =
            SymbolicTensor::contract_ports(&left, &right, &[PortPair { left: 0, right: 0 }])
                .unwrap();
        let indexed = SymbolicTensor::materialize_interface_ports(
            &contracted,
            &HashMap::from([(0, AbstractIndex::Normal(23))]),
        )
        .unwrap();

        let left_slots = factor_slots(&indexed, left_symbol);
        let right_slots = factor_slots(&indexed, right_symbol);
        assert_eq!(left_slots.len(), 2);
        assert_eq!(right_slots.len(), 1);
        assert_eq!(left_slots[0].aind(), right_slots[0].aind());
        assert_ne!(left_slots[0].aind(), AbstractIndex::Normal(23));
        assert_eq!(left_slots[1].aind(), AbstractIndex::Normal(23));
    }

    #[test]
    fn open_contractions_allocate_fresh_dummies() {
        let ports = [
            rep().slot(PartialIndex::open(0)),
            rep().slot(PartialIndex::open(1)),
        ];
        let left_symbol = SPENSO_TAG.tensor_symbol("fresh_dummy_left");
        let right_symbol = SPENSO_TAG.tensor_symbol("fresh_dummy_right");
        let left = partial_tensor(left_symbol, &ports, &ports);
        let right = partial_tensor(right_symbol, &ports, &ports);
        let first =
            SymbolicTensor::contract_ports(&left, &right, &[PortPair { left: 0, right: 0 }])
                .unwrap();
        let second =
            SymbolicTensor::contract_ports(&left, &right, &[PortPair { left: 0, right: 0 }])
                .unwrap();

        let first_dummy = factor_slots(&first.expression, left_symbol)[0].aind();
        let first_partner = factor_slots(&first.expression, right_symbol)[0].aind();
        let second_dummy = factor_slots(&second.expression, left_symbol)[0].aind();
        let second_partner = factor_slots(&second.expression, right_symbol)[0].aind();
        assert_eq!(first_dummy, first_partner);
        assert_eq!(second_dummy, second_partner);
        assert_ne!(first_dummy, second_dummy);
    }

    #[test]
    fn open_contraction_dummy_does_not_alias_an_explicit_dummy_symbol() {
        let next_dummy = match AbstractIndex::new_dummy() {
            AbstractIndex::Dummy(index) => AbstractIndex::Dummy(index + 1),
            _ => unreachable!(),
        };
        let explicit = AbstractIndex::try_from(next_dummy.to_atom().as_view()).unwrap();
        let open = rep().slot(PartialIndex::open(0));
        let spectator = rep().slot(PartialIndex::Explicit(explicit));
        let left_symbol = SPENSO_TAG.tensor_symbol("dummy_collision_left");
        let right_symbol = SPENSO_TAG.tensor_symbol("dummy_collision_right");
        let left = partial_tensor(left_symbol, &[open, spectator], &[open, spectator]);
        let right = partial_tensor(right_symbol, &[open], &[open]);

        let contracted =
            SymbolicTensor::contract_ports(&left, &right, &[PortPair { left: 0, right: 0 }])
                .unwrap();
        let left_slots = factor_slots(&contracted.expression, left_symbol);
        let right_slots = factor_slots(&contracted.expression, right_symbol);

        assert_eq!(left_slots[0].aind(), right_slots[0].aind());
        assert_ne!(left_slots[0].aind(), explicit);
        assert_eq!(left_slots[1].aind(), explicit);
    }

    #[test]
    fn explicit_indices_must_agree_before_contraction() {
        let left_port = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(3)));
        let right_port = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(4)));
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("unequal_explicit_left"),
            &[left_port],
            &[left_port],
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("unequal_explicit_right"),
            &[right_port],
            &[right_port],
        );

        assert!(matches!(
            SymbolicTensor::contract_ports(&left, &right, &[PortPair { left: 0, right: 0 }]),
            Err(TensorCompositionError::UnequalExplicitIndices { left: 0, right: 0 })
        ));
    }

    #[test]
    fn index_materialization_rejects_a_third_explicit_occurrence() {
        let representation = rep();
        let index = AbstractIndex::Normal(27);
        let open = representation.slot(PartialIndex::open(0));
        let explicit = representation.slot(PartialIndex::Explicit(index));
        let value = partial_tensor(
            SPENSO_TAG.tensor_symbol("third_index_occurrence"),
            &[open, explicit, explicit],
            &[open],
        );

        assert!(matches!(
            SymbolicTensor::materialize_interface_ports(&value, &HashMap::from([(0, index)])),
            Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: actual,
                occurrences: 3
            }) if actual == index
        ));
    }

    #[test]
    fn explicit_contract_rejects_a_third_explicit_occurrence() {
        let representation = rep();
        let index = AbstractIndex::Normal(31);
        let open = representation.slot(PartialIndex::open(0));
        let explicit = representation.slot(PartialIndex::Explicit(index));
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("third_contract_occurrence_left"),
            &[open, explicit],
            &[open, explicit],
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("third_contract_occurrence_right"),
            &[explicit],
            &[explicit],
        );

        assert!(matches!(
            SymbolicTensor::contract_ports(&left, &right, &[PortPair { left: 0, right: 0 }]),
            Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: actual,
                occurrences: 3
            }) if actual == index
        ));
    }

    #[test]
    fn explicit_compose_rejects_a_third_explicit_occurrence() {
        let representation = rep();
        let index = AbstractIndex::Normal(33);
        let explicit = representation.slot(PartialIndex::Explicit(index));
        let left_ports = [
            representation.slot(PartialIndex::open(0)),
            representation.slot(PartialIndex::open(1)),
            explicit,
        ];
        let right_ports = [explicit, representation.slot(PartialIndex::open(0))];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("third_compose_occurrence_left"),
            &left_ports,
            &left_ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("third_compose_occurrence_right"),
            &right_ports,
            &right_ports,
        );

        assert!(matches!(
            SymbolicTensor::compose(
                &left,
                &right,
                MatrixChannel {
                    input: 0,
                    output: 1,
                },
                MatrixChannel {
                    input: 0,
                    output: 1,
                },
            ),
            Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: actual,
                occurrences: 3
            }) if actual == index
        ));
    }

    #[test]
    fn explicit_occurrences_are_counted_per_additive_branch() {
        let representation = rep();
        let index = AbstractIndex::Normal(35);
        let explicit = representation.slot(PartialIndex::Explicit(index));
        let terms = ["first", "second", "third"].map(|name| {
            partial_tensor(SPENSO_TAG.tensor_symbol(name), &[explicit], &[explicit]).expression
        });
        let sum = terms.into_iter().fold(Atom::Zero, |sum, term| sum + term);

        validate_explicit_index_occurrences(&sum).unwrap();
    }

    #[test]
    fn occurrence_summaries_reuse_callback_leaves_without_replaying_them() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };

        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = tensor_symbol!(
            "occurrence_summary_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let slot = rep().slot::<AbstractIndex, _>(AbstractIndex::Normal(73));
        let leaf = FunctionBuilder::new(head).add_arg(slot.to_atom()).finish();
        let [a, b, c, d] = ["summary_a", "summary_b", "summary_c", "summary_d"].map(|name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(&leaf)
                .finish()
        });
        let valid = (a + b) * c;
        let invalid = &valid * d;
        calls.store(0, Ordering::Relaxed);

        let mut slots = SlotMatcher::default();
        let mut cache = AHashMap::new();
        let mut occurrences = ExplicitIndexOccurrences::default();
        OCCURRENCE_FUNCTION_VISITS.set(0);
        occurrences.append(&InterfaceInference::index_occurrences(
            leaf.as_view(),
            &mut slots,
            &mut cache,
        ));
        let visits = OCCURRENCE_FUNCTION_VISITS.get();
        assert!(visits > 0);
        for expected in [2, 3] {
            occurrences.append(&InterfaceInference::index_occurrences(
                leaf.as_view(),
                &mut slots,
                &mut cache,
            ));
            assert_eq!(occurrences.count(slot), expected);
            assert_eq!(OCCURRENCE_FUNCTION_VISITS.get(), visits);
        }
        validate_explicit_index_occurrences(&valid).unwrap();
        assert!(matches!(
            validate_explicit_index_occurrences(&invalid),
            Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: AbstractIndex::Normal(73),
                occurrences: 3
            })
        ));
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn occurrence_summaries_cache_scalar_dots_and_keep_encoded_index_identity() {
        let compact = rep().to_symbolic([]);
        let p = FunctionBuilder::new(vector_symbol!(occurrence_summary_p))
            .add_arg(&compact)
            .finish();
        let q = FunctionBuilder::new(vector_symbol!(occurrence_summary_q))
            .add_arg(&compact)
            .finish();
        let dot = FunctionBuilder::new(ETS.metric).add_args([p, q]).finish();
        let mut slots = SlotMatcher::default();
        let mut cache = AHashMap::new();
        let mut occurrences = ExplicitIndexOccurrences::default();
        OCCURRENCE_FUNCTION_VISITS.set(0);
        occurrences.append(&InterfaceInference::index_occurrences(
            dot.as_view(),
            &mut slots,
            &mut cache,
        ));
        let visits = OCCURRENCE_FUNCTION_VISITS.get();
        assert!(visits > 0);
        assert!(occurrences.is_empty());
        occurrences.append(&InterfaceInference::index_occurrences(
            dot.as_view(),
            &mut slots,
            &mut cache,
        ));
        assert_eq!(OCCURRENCE_FUNCTION_VISITS.get(), visits);
        assert!(occurrences.is_empty());

        let head = SPENSO_TAG.tensor_symbol("occurrence_summary_open_owner");
        let first = rep().slot::<AbstractIndex, _>(AbstractIndex::Open {
            owner: 101,
            axis: 0,
        });
        let second = rep().slot::<AbstractIndex, _>(AbstractIndex::Open {
            owner: 102,
            axis: 0,
        });
        for slot in [first, second, first] {
            let leaf = FunctionBuilder::new(head).add_arg(slot.to_atom()).finish();
            occurrences.append(&InterfaceInference::index_occurrences(
                leaf.as_view(),
                &mut slots,
                &mut cache,
            ));
        }
        assert_eq!(occurrences.count(first), 2);
        assert_eq!(occurrences.count(second), 1);
    }

    #[test]
    fn explicit_multiplicity_preserves_duality_dimensions_and_metadata() {
        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let index = AbstractIndex::Normal(39);
        let base_slot = base.slot::<AbstractIndex, _>(index).to_atom();
        let dual_slot = base.dual().slot::<AbstractIndex, _>(index).to_atom();
        let metadata = FunctionBuilder::new(symbolica::symbol!("multiplicity_metadata"))
            .add_args([&base_slot, &dual_slot, &base_slot])
            .finish();
        let pair = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("multiplicity_pair"))
            .add_arg(&metadata)
            .add_arg(&base_slot)
            .add_arg(&dual_slot)
            .finish();
        let other_dimension: Representation<LibraryRep> = ColorFundamental {}.new_rep(4).cast();
        let unrelated = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("multiplicity_other_dim"))
            .add_arg(other_dimension.slot::<AbstractIndex, _>(index).to_atom())
            .finish();
        let valid = pair * metadata * unrelated;
        validate_explicit_index_occurrences(&valid).unwrap();

        let third = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("multiplicity_third"))
            .add_arg(base_slot)
            .finish();
        assert!(matches!(
            validate_explicit_index_occurrences(&(valid * third)),
            Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: actual,
                occurrences: 3
            }) if actual == index
        ));
    }

    #[test]
    fn explicit_multiplicity_keeps_nested_sum_and_power_scopes() {
        let index = AbstractIndex::Normal(41);
        let slot = rep().slot(PartialIndex::Explicit(index));
        let [a, b, c, d] = ["scope_a", "scope_b", "scope_c", "scope_d"].map(|name| {
            partial_tensor(SPENSO_TAG.tensor_symbol(name), &[slot], &[slot]).expression
        });
        let alternatives = (&a + &b) * &c;
        validate_explicit_index_occurrences(&alternatives).unwrap();
        assert!(matches!(
            validate_explicit_index_occurrences(&(alternatives * d)),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));

        // Count literal open ports without materializing repeated factors.
        assert!(matches!(
            validate_explicit_index_occurrences(&(a.pow(3) * &b)),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { .. })
        ));
        validate_explicit_index_occurrences(&(a.pow(&b) * &c)).unwrap();
        let exponent = a * b * c;
        let power = Atom::var(symbolica::symbol!("multiplicity_scalar_base")).pow(exponent);
        assert!(matches!(
            validate_explicit_index_occurrences(&power),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));
    }

    #[test]
    fn explicit_multiplicity_keeps_different_branch_indices_separate() {
        let ports =
            [43, 47].map(|index| rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(index))));
        let [a, b, c, d] = [
            ("branch_a", ports[0]),
            ("branch_b", ports[1]),
            ("branch_c", ports[0]),
            ("branch_d", ports[1]),
        ]
        .map(|(name, port)| {
            partial_tensor(SPENSO_TAG.tensor_symbol(name), &[port], &[port]).expression
        });
        let alternatives = (&a + b) * (c + d);
        validate_explicit_index_occurrences(&alternatives).unwrap();

        // A third occurrence in only one possible branch still invalidates it.
        let extra = partial_tensor(
            SPENSO_TAG.tensor_symbol("branch_extra"),
            &ports[..1],
            &ports[..1],
        );
        assert!(matches!(
            validate_explicit_index_occurrences(&(alternatives * extra.expression)),
            Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: AbstractIndex::Normal(43),
                occurrences: 3
            })
        ));
    }

    #[test]
    fn materialization_reuses_the_selected_index_in_every_summand() {
        let port = rep().slot(PartialIndex::open(0));
        let first = partial_tensor(
            SPENSO_TAG.tensor_symbol("sum_materialization_first"),
            &[port],
            &[port],
        );
        let second = partial_tensor(
            SPENSO_TAG.tensor_symbol("sum_materialization_second"),
            &[port],
            &[port],
        );
        let sum = SymbolicTensor::new(
            first.expression + second.expression,
            PartialStructure::from_logical_slots([port]),
        );
        let materialized = SymbolicTensor::materialize_interface_ports(
            &sum,
            &HashMap::from([(0, AbstractIndex::Normal(29))]),
        )
        .unwrap();

        let AtomView::Add(sum) = materialized.as_view() else {
            panic!("expected a sum")
        };
        assert_eq!(sum.get_nargs(), 2);
        for term in sum.iter() {
            let AtomView::Fun(tensor) = term else {
                panic!("expected a tensor summand")
            };
            let slot =
                Slot::<LibraryRep, AbstractIndex>::try_from(tensor.iter().next().unwrap()).unwrap();
            assert_eq!(slot.aind(), AbstractIndex::Normal(29));
        }
    }

    #[test]
    fn matrix_roles_follow_self_dual_order_and_dual_orientation() {
        let self_dual = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let self_dual_ports = [
            spectator.slot(PartialIndex::open(0)),
            self_dual.slot(PartialIndex::open(1)),
            self_dual.slot(PartialIndex::open(2)),
        ];
        let self_dual_tensor = partial_tensor(
            SPENSO_TAG.tensor_symbol("self_dual_channel"),
            &self_dual_ports,
            &self_dual_ports,
        );
        assert_eq!(
            SymbolicTensor::matrix_channel(&self_dual_tensor),
            Some(MatrixChannel {
                input: 1,
                output: 2
            })
        );

        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let dual = base.dual();
        let dual_ports = [
            dual.slot(PartialIndex::open(0)),
            spectator.slot(PartialIndex::open(1)),
            base.slot(PartialIndex::open(2)),
        ];
        let dual_tensor = partial_tensor(
            SPENSO_TAG.tensor_symbol("dual_channel"),
            &dual_ports,
            &dual_ports,
        );
        assert_eq!(
            SymbolicTensor::matrix_channel(&dual_tensor),
            Some(MatrixChannel {
                input: 2,
                output: 0
            })
        );
    }

    #[test]
    fn gamma_and_color_chain_factors_keep_registered_argument_order() {
        let mink: Representation<LibraryRep> = Minkowski {}.new_rep(4).cast();
        let bis: Representation<LibraryRep> = Bispinor {}.new_rep(4).cast();
        let gamma_logical = [
            mink.slot(PartialIndex::open(0)),
            bis.slot(PartialIndex::open(1)),
            bis.slot(PartialIndex::open(2)),
        ];
        let gamma = partial_tensor(
            AGS.gamma,
            &[gamma_logical[1], gamma_logical[2], gamma_logical[0]],
            &gamma_logical,
        );
        let gamma_factor =
            SymbolicTensor::chain_factors(&gamma, SymbolicTensor::matrix_channel(&gamma).unwrap())
                .unwrap()
                .pop()
                .unwrap();
        let AtomView::Fun(gamma_factor) = gamma_factor.as_view() else {
            panic!("expected a gamma factor")
        };
        let gamma_args = gamma_factor.iter().collect::<Vec<_>>();
        assert!(matches!(gamma_args[0], AtomView::Var(v) if v.get_symbol() == SPENSO_TAG.chain_in));
        assert!(
            matches!(gamma_args[1], AtomView::Var(v) if v.get_symbol() == SPENSO_TAG.chain_out)
        );
        assert_eq!(
            Representation::<LibraryRep>::try_from(gamma_args[2]).unwrap(),
            mink
        );

        let adjoint: Representation<LibraryRep> = ColorAdjoint {}.new_rep(8).cast();
        let fundamental: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let antifundamental = fundamental.dual();
        let color_ports = [
            adjoint.slot(PartialIndex::open(0)),
            fundamental.slot(PartialIndex::open(1)),
            antifundamental.slot(PartialIndex::open(2)),
        ];
        let color = partial_tensor(CS.t, &color_ports, &color_ports);
        let color_factor =
            SymbolicTensor::chain_factors(&color, SymbolicTensor::matrix_channel(&color).unwrap())
                .unwrap()
                .pop()
                .unwrap();
        let AtomView::Fun(color_factor) = color_factor.as_view() else {
            panic!("expected a color-generator factor")
        };
        let color_args = color_factor.iter().collect::<Vec<_>>();
        assert_eq!(
            Representation::<LibraryRep>::try_from(color_args[0]).unwrap(),
            adjoint
        );
        assert!(matches!(color_args[1], AtomView::Var(v) if v.get_symbol() == SPENSO_TAG.chain_in));
        assert!(
            matches!(color_args[2], AtomView::Var(v) if v.get_symbol() == SPENSO_TAG.chain_out)
        );
    }

    #[test]
    fn tracing_a_chain_spectator_keeps_its_original_endpoints() {
        let channel_rep = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let ports = [
            channel_rep.slot(PartialIndex::open(0)),
            channel_rep.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::open(2)),
        ];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("spectator_trace_left"),
            &ports,
            &ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("spectator_trace_right"),
            &ports,
            &ports,
        );
        let chain = SymbolicTensor::compose(
            &left,
            &right,
            MatrixChannel {
                input: 0,
                output: 1,
            },
            MatrixChannel {
                input: 0,
                output: 1,
            },
        )
        .unwrap();
        let traced = SymbolicTensor::trace_ports(
            &chain,
            MatrixChannel {
                input: 2,
                output: 3,
            },
        )
        .unwrap();

        assert_eq!(traced.rank(), 2);
        assert!(
            matches!(traced.expression.as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.chain)
        );
        traced.expression.validate_chain_like_nesting().unwrap();
        let materialized = SymbolicTensor::materialize_interface_ports(
            &traced,
            &HashMap::from([
                (0, AbstractIndex::Normal(31)),
                (1, AbstractIndex::Normal(37)),
            ]),
        )
        .unwrap();
        let network = materialized
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert_eq!(network.graph.dangling_indices().len(), 2);
    }

    #[test]
    fn tracing_chain_endpoints_then_spectators_keeps_one_root_trace() {
        let channel_rep = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let ports = [
            channel_rep.slot(PartialIndex::open(0)),
            channel_rep.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::open(2)),
        ];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("double_trace_left"),
            &ports,
            &ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("double_trace_right"),
            &ports,
            &ports,
        );
        let chain = SymbolicTensor::compose(
            &left,
            &right,
            MatrixChannel {
                input: 0,
                output: 1,
            },
            MatrixChannel {
                input: 0,
                output: 1,
            },
        )
        .unwrap();

        let primary = SymbolicTensor::trace_ports(
            &chain,
            MatrixChannel {
                input: 0,
                output: 1,
            },
        )
        .unwrap();
        let traced = SymbolicTensor::trace_unique(&primary).unwrap();

        assert!(traced.is_scalar());
        assert!(
            matches!(traced.expression.as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.trace)
        );
        traced.expression.validate_chain_like_nesting().unwrap();
        let network = traced
            .expression
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(network.state.is_scalar());
        assert!(network.graph.dangling_indices().is_empty());
    }

    #[test]
    fn trace_interface_follows_canonical_cyclic_factor_rotation() {
        let channel_rep: Representation<LibraryRep> = Bispinor {}.new_rep(4).cast();
        let spectator: Representation<LibraryRep> = Minkowski {}.new_rep(4).cast();
        let left_index = AbstractIndex::Normal(43);
        let right_index = AbstractIndex::Normal(41);
        let left_ports = [
            spectator.slot(PartialIndex::Explicit(left_index)),
            channel_rep.slot(PartialIndex::open(0)),
            channel_rep.slot(PartialIndex::open(1)),
        ];
        let right_ports = [
            spectator.slot(PartialIndex::Explicit(right_index)),
            channel_rep.slot(PartialIndex::open(0)),
            channel_rep.slot(PartialIndex::open(1)),
        ];
        let left = partial_tensor(
            AGS.gamma,
            &[left_ports[1], left_ports[2], left_ports[0]],
            &left_ports,
        );
        let right = partial_tensor(
            AGS.gamma,
            &[right_ports[1], right_ports[2], right_ports[0]],
            &right_ports,
        );

        let traced = SymbolicTensor::trace_unique(
            &SymbolicTensor::compose(
                &left,
                &right,
                MatrixChannel {
                    input: 1,
                    output: 2,
                },
                MatrixChannel {
                    input: 1,
                    output: 2,
                },
            )
            .unwrap(),
        )
        .unwrap();

        assert_eq!(
            traced
                .structure
                .logical_slots()
                .into_iter()
                .map(|slot| slot.aind)
                .collect::<Vec<_>>(),
            vec![
                PartialIndex::Explicit(right_index),
                PartialIndex::Explicit(left_index),
            ]
        );
    }

    #[test]
    fn trace_closes_chain_in_canonical_form() {
        let first = tensor("first", &[rep(), rep()]);
        let second = tensor("second", &[rep(), rep()]);
        let chain = first.multiply(&second).unwrap();
        let traced = SymbolicTensor::trace_unique(&chain).unwrap();

        assert!(traced.is_scalar());
        assert!(
            matches!(traced.expression.as_view(), AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.trace)
        );
        assert!(
            shadowing::trace_parts(match traced.expression.as_view() {
                AtomView::Fun(fun) => fun,
                _ => unreachable!(),
            })
            .is_some()
        );
    }

    #[test]
    fn explicit_channels_require_distinct_endpoints() {
        let matrix = tensor("degenerate_channel", &[rep(), rep()]);
        let channel = MatrixChannel {
            input: 0,
            output: 0,
        };

        assert!(matches!(
            SymbolicTensor::trace_ports(&matrix, channel),
            Err(TensorCompositionError::DegenerateChannel { .. })
        ));
        assert!(matches!(
            SymbolicTensor::compose(
                &matrix,
                &matrix,
                channel,
                MatrixChannel {
                    input: 0,
                    output: 1,
                },
            ),
            Err(TensorCompositionError::DegenerateChannel { .. })
        ));
    }

    #[test]
    fn explicit_channels_validate_representation_and_dual_orientation() {
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let incompatible_left = tensor("incompatible_channel_left", &[rep(), minkowski]);
        let incompatible_right = tensor("incompatible_channel_right", &[minkowski, rep()]);
        let channel = MatrixChannel {
            input: 0,
            output: 1,
        };
        assert!(matches!(
            SymbolicTensor::compose(&incompatible_left, &incompatible_right, channel, channel),
            Err(TensorCompositionError::IncompatiblePorts { .. })
        ));

        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let reversed = tensor("reversed_dual_channel", &[base.dual(), base]);
        assert!(matches!(
            SymbolicTensor::trace_ports(&reversed, channel),
            Err(TensorCompositionError::InvalidChannelOrientation {
                input: 0,
                output: 1
            })
        ));
    }
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(
        expected = "normalization callback contract violated in reindex_interface_ports"
    )]
    fn relabeling_checks_callback_contract() {
        let zero = Arc::new(AtomicBool::new(false));
        let emit_zero = Arc::clone(&zero);
        let target = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(419)));
        let target_atom = port_atom(target);
        let head = tensor_symbol!(
            "relabel_result_rank_violation",
            norm = move |function, output| {
                if let AtomView::Fun(function) = function
                    && function.iter().any(|arg| arg == target_atom.as_view())
                {
                    **output = if emit_zero.load(Ordering::Relaxed) {
                        Atom::Zero
                    } else {
                        Atom::num(7)
                    };
                }
            }
        );
        let port = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(421)));
        let value = partial_tensor(head, &[port], &[port]);
        assert!(!PortRewriteState::exact_arithmetic(
            value.expression.as_view()
        ));
        let replacements = HashMap::from([(0, AbstractIndex::Normal(419))]);
        // The relabelling operation determines the target interface. Normalizer
        // output does not replace admitted tensor metadata, even when scalar.
        value.reindex_interface_ports(&replacements).unwrap();
    }

    #[test]
    fn positional_ranked_zero_seeds_declared_index_history() {
        let i = AbstractIndex::Normal(99801);
        let j = AbstractIndex::Normal(99803);
        let ports = [
            rep().slot(PartialIndex::Explicit(i)),
            rep().slot(PartialIndex::open(1)),
        ];
        let vector = vector("positional_zero_vector", ports[0]);
        for admitted in [false, true] {
            let zero = SymbolicTensor::new(Atom::Zero, PartialStructure::from_logical_slots(ports));
            if admitted {
                zero.ensure_validated().unwrap();
            }
            let contracted = zero
                .contract_ports(&vector, &[PortPair { left: 0, right: 0 }])
                .unwrap();
            assert!(contracted.is_zero());
            assert_eq!(contracted.rank(), 1);
            let reindexed = zero
                .reindex_interface_ports(&HashMap::from([(0, j)]))
                .unwrap();
            assert!(reindexed.is_zero());
            assert_eq!(reindexed.rank(), 2);
            assert_eq!(
                reindexed
                    .explicit_index_occurrences()
                    .unwrap()
                    .count(rep().slot(j)),
                1
            );
        }
    }

    #[test]
    fn multiply_accepts_raw_positional_logical_order() {
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let physical = [
            rep().slot(PartialIndex::open(0)),
            mink.slot(PartialIndex::open(1)),
        ];
        let logical = [
            mink.slot(PartialIndex::open(0)),
            rep().slot(PartialIndex::open(1)),
        ];
        let source = partial_tensor(
            tensor_symbol!("raw_positional_product"),
            &physical,
            &logical,
        );
        let vector = tensor("raw_positional_vector", &[mink]);
        let direct = source
            .contract_ports(&vector, &[PortPair { left: 0, right: 0 }])
            .unwrap();
        let automatic = source.multiply(&vector).unwrap();
        assert_eq!(automatic.structure, direct.structure);
        assert_eq!(automatic.rank(), 1);
        assert!((&source * &vector).is_ok());
        assert_eq!(
            source.try_scalar_mul(&Atom::num(3)).unwrap().structure,
            source.structure
        );
    }

    #[test]
    fn scalar_multiplication_keeps_closed_dummy_scope_independent() {
        let index = AbstractIndex::Normal(99811);
        let port = rep().slot(PartialIndex::Explicit(index));
        let source = partial_tensor(tensor_symbol!("scalar_scope_vector"), &[port], &[port]);
        let scalar = function!(tensor_symbol!("scalar_scope_p"), port_atom(port))
            * function!(tensor_symbol!("scalar_scope_q"), port_atom(port));
        let result = source
            .scalar_mul(&scalar)
            .expect("a closed scalar must multiply a tensor");
        assert_eq!(result.structure, source.structure);
        SymbolicTensor::validate_interface(&result.expression, &result.structure).unwrap();
        let indexed = source.clone().try_into_indexed().unwrap();
        let scaled = indexed
            .scalar_mul(&scalar)
            .expect("indexed scalar multiplication must succeed");
        assert_eq!(scaled.structure(), indexed.structure());
    }
}
