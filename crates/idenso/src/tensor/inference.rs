//! Tensor-aware interface inference shared by Rust and Python entry points.
//!
//! Cached interfaces retain raw logical ports. Only enclosing algebra contracts
//! explicit pairs; unresolved port identities are canonicalized per occurrence.
//! See the [admission example](SymbolicTensor).
use super::composition::{ExplicitIndexOccurrences, TensorCompositionError};
use super::raw_leaves::RawLeafInterfaces;
use super::regions::AdmittedInterfaces;
use super::simplification::observation::DomainObservations;
use super::{SymbolicTensor, TensorProofs, composition};
use crate::color::CS;
use crate::dirac::AGS;
use crate::epsilon::EPSILON_SYMBOL;
#[cfg(feature = "reference-cases")]
use crate::reference_cases::timing;
#[cfg(feature = "reference-cases")]
use crate::reference_cases::timing::Phase;
use crate::representations::{Bispinor, ColorAdjoint, ColorAntiFundamental, ColorFundamental};
use crate::shorthands::schoonschip::{DotNormalizer, SlotContraction};
use crate::{CookMode, CookSettings, CookingError, NetworkToolingError};
use ahash::{AHashMap, AHashSet};
use spenso::network::TensorNetworkError;
use spenso::network::library::DummyKey;
use spenso::network::library::symbolic::{ETS, ExplicitKey};
use spenso::network::parsing::{
    AtomStructureExt, LeafInterfaceCache, ParseState, StrictTensorFilter, StructureFromAtom,
};
use spenso::network::tags::SPENSO_TAG;
use spenso::shadowing;
use spenso::shadowing::{ANTISYM, CYCLIC, SYM};
#[cfg(debug_assertions)]
use spenso::structure::StructureContract;
use spenso::structure::abstract_index::{AIND_SYMBOLS, AbstractIndex};
use spenso::structure::dimension::Dimension;
use spenso::structure::partial::{
    PartialIndex, PartialSlot, PartialStructure, PartialStructureExt,
};
use spenso::structure::representation::{LibraryRep, Minkowski, RepName, Representation};
use spenso::structure::slot::{
    AbsInd, IsAbstractSlot, ParseableAind, Slot, SlotError, SlotMatch, SlotMatcher, SlotView,
};
use spenso::structure::{
    Canonicalized, IndexedStructure, OrderedStructure, StructureError, TensorStructure,
};
#[cfg(all(test, debug_assertions))]
use std::cell::Cell;
use std::cell::{OnceCell, RefCell};
use std::collections::{HashMap, HashSet};
use std::mem;
use std::slice;
use std::sync::{Arc, Mutex, OnceLock};
use symbolica::atom::{
    Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, NamespacedSymbol, Symbol,
    SymbolAttribute, SymbolBuilder,
};
use symbolica::domains::rational::Rational;
use symbolica::evaluate::EvaluatorBuilder;
use thiserror::Error;

/// Failures of tensor admission, metadata transport, and typed arithmetic.
///
/// See the [admission example](SymbolicTensor).
#[derive(Debug, Error)]
pub enum TensorInferenceError {
    #[error("transformation interrupted")]
    Interrupted,
    #[error("{0}")]
    Invalid(String),
    #[error("predefined tensor `{factory}` requires {signature}")]
    InvalidBuiltinSignature {
        factory: &'static str,
        signature: &'static str,
    },
    #[error(transparent)]
    Composition(#[from] TensorCompositionError),
    #[error(transparent)]
    Network(#[from] NetworkToolingError),
    /// The expression could not be represented by the planning graph.
    #[error(transparent)]
    Planning(Box<TensorNetworkError<DummyKey, Symbol>>),
    #[error(transparent)]
    Cooking(#[from] CookingError),
    #[error("cannot encode mapped dimension: {0}")]
    MappedDimensionEncoding(#[source] CookingError),
    #[error("{operation}: {source}")]
    Slot {
        operation: &'static str,
        #[source]
        source: SlotError,
    },
}

impl From<TensorNetworkError<DummyKey, Symbol>> for TensorInferenceError {
    fn from(source: TensorNetworkError<DummyKey, Symbol>) -> Self {
        Self::Planning(Box::new(source))
    }
}

impl TensorInferenceError {
    pub(crate) fn invalid(message: impl Into<String>) -> Self {
        Self::Invalid(message.into())
    }
}

type InferenceResult<T> = Result<T, TensorInferenceError>;

mod observation;
pub(super) use observation::ObservationScope;

/// A transient view of one admitted source. Its proof cells are borrowed so
/// observations update the original carrier rather than a detached clone.
#[derive(Clone, Copy)]
pub(crate) struct BorrowedTensorParts<'a> {
    pub(super) expression: &'a Atom,
    pub(super) structure: &'a PartialStructure,
    pub(super) is_metric: bool,
    pub(super) is_composite: bool,
    pub(super) proofs: &'a TensorProofs,
    /// Original input declarations, independent of a materializer's generated facts.
    pub(super) source_proofs: &'a TensorProofs,
    pub(super) work: Option<&'a super::simplification::ReductionCache>,
}

/// Project stored wrapping interfaces without rebuilding index atoms. Generic
/// index types are matched against existing result slots; metadata remains the
/// authority for their identities, including the nominal slots of zero.
pub(super) trait TensorInterface {
    /// Transfer the current free interface to algebra, retaining positional
    /// bindings while dropping a named leaf's admission-layout history.
    fn into_algebra_result(self) -> Self;

    #[cfg(debug_assertions)]
    fn to_partial_interface(&self, expression: &Atom) -> InferenceResult<PartialStructure>;
}

impl TensorInterface for PartialStructure {
    fn into_algebra_result(self) -> Self {
        if SymbolicTensor::interface_is_named(&self) {
            Canonicalized::identity(self.into_canonical())
        } else {
            self
        }
    }

    #[cfg(debug_assertions)]
    fn to_partial_interface(&self, _expression: &Atom) -> InferenceResult<PartialStructure> {
        Ok(self.clone())
    }
}

impl TensorInterface for IndexedStructure {
    fn into_algebra_result(self) -> Self {
        if self
            .canonical()
            .external_indices_iter()
            .all(|index| SymbolicTensor::index_is_named(PartialIndex::Explicit(index)))
        {
            Canonicalized::identity(self.into_canonical())
        } else {
            self
        }
    }

    #[cfg(debug_assertions)]
    fn to_partial_interface(&self, _expression: &Atom) -> InferenceResult<PartialStructure> {
        Ok(self.clone().into_partial())
    }
}

impl<A: AbsInd + ParseableAind> TensorInterface for OrderedStructure<LibraryRep, A> {
    fn into_algebra_result(self) -> Self {
        self
    }

    #[cfg(debug_assertions)]
    fn to_partial_interface(&self, expression: &Atom) -> InferenceResult<PartialStructure> {
        let mut interface = self.clone();
        interface.trace_out();
        if expression.is_zero() {
            // Zero has no syntax witnessing labels; its stored shape suffices
            // because the contract permits zero with every nominal interface.
            return Ok(PartialStructure::from_logical_slots(
                interface
                    .external_reps_iter()
                    .enumerate()
                    .map(|(axis, rep)| rep.slot(PartialIndex::open(axis))),
            ));
        }
        let mut indices = interface
            .external_structure_iter()
            .map(|slot| (slot, None::<AbstractIndex>))
            .collect::<HashMap<_, _>>();
        expression.visitor(&mut |node| {
            if let Ok(slot) = Slot::<LibraryRep, A>::try_from(node)
                && let Some(index) = indices.get_mut(&slot)
                && let Ok(observed) = Slot::<LibraryRep, AbstractIndex>::try_from(node)
            {
                *index = Some(observed.aind());
            }
            true
        });
        let slots = interface
            .external_structure_iter()
            .map(|slot| {
                let index = indices[&slot].ok_or_else(|| {
                    TensorInferenceError::invalid(
                        "wrapped result does not contain an admitted free index",
                    )
                })?;
                Ok(slot.rep().slot(PartialIndex::Explicit(index)))
            })
            .collect::<InferenceResult<Vec<_>>>()?;
        Ok(PartialStructure::from_logical_slots(slots))
    }
}

/// Separate debug verification from production work in inference budgets.
#[cfg(all(test, debug_assertions))]
#[derive(Debug)]
pub(crate) struct DebugInference {
    previous: bool,
}

#[cfg(all(test, debug_assertions))]
thread_local! {
    static DEBUG_INFERENCE_ACTIVE: Cell<bool> = const { Cell::new(false) };
}

#[cfg(all(test, debug_assertions))]
impl DebugInference {
    pub(crate) fn enter() -> Self {
        Self {
            previous: DEBUG_INFERENCE_ACTIVE.with(|active| active.replace(true)),
        }
    }

    pub(crate) fn is_active() -> bool {
        DEBUG_INFERENCE_ACTIVE.with(Cell::get)
    }
}

#[cfg(all(test, debug_assertions))]
impl Drop for DebugInference {
    fn drop(&mut self) {
        DEBUG_INFERENCE_ACTIVE.with(|active| active.set(self.previous));
    }
}

#[derive(Debug)]
struct ObservedTensor {
    interface: PartialStructure,
    occurrences: Option<ExplicitIndexOccurrences>,
    interfaces: AdmittedInterfaces,
    leaf_interfaces: LeafInterfaceCache,
}

impl TensorProofs {
    pub(crate) fn owns_published_interface(
        &self,
        expression: AtomView<'_>,
        interface: &PartialStructure,
    ) -> bool {
        self.owns_expression(expression) && self.published_interface_matches(expression, interface)
    }

    /// Check the user-normalizer contract without replaying normalization.
    /// Release builds trust the producer; zero keeps its nominal interface.
    /// `rebuilt_head` identifies a newly applied head even if its callback removes it.
    /// Explicit encoded indices retain their owner/axis identities during observation.
    pub(crate) fn assert_normalization_preserves_interface(
        operation: &str,
        sources: &[&Atom],
        expression: &Atom,
        structure: Result<&PartialStructure, &TensorInferenceError>,
        rebuilt_head: Option<Symbol>,
    ) {
        #[cfg(all(test, debug_assertions))]
        let _debug = DebugInference::enter();
        #[cfg(debug_assertions)]
        if rebuilt_head.is_some_and(|head| !InterfaceInference::intrinsic_normalization_head(head))
            || sources
                .iter()
                .any(|source| !InterfaceInference::normalization_is_intrinsic(source.as_view()))
        {
            let source = sources
                .iter()
                .find(|source| !InterfaceInference::normalization_is_intrinsic(source.as_view()))
                .or_else(|| sources.first());
            Self::assert_interface_agreement(
                operation,
                source.copied(),
                rebuilt_head,
                expression,
                structure,
            );
        }
        #[cfg(not(debug_assertions))]
        let _ = (operation, sources, expression, structure, rebuilt_head);
    }
    #[cfg(debug_assertions)]
    fn assert_interface_agreement(
        operation: &str,
        source: Option<&Atom>,
        rebuilt_head: Option<Symbol>,
        expression: &Atom,
        structure: Result<&PartialStructure, &TensorInferenceError>,
    ) {
        let checked = structure.ok().map(|structure| {
            InterfaceInference::merge_explicit_interface_sequence(slice::from_ref(structure))
                .and_then(|external| {
                    // Opaque scalar admission establishes an empty interface without
                    // validating tensor-free payload syntax. Check the same contract.
                    if external.canonical().is_scalar()
                        && !InterfaceInference::contains_tensor_symbol(expression.as_view())
                    {
                        return Ok(());
                    }
                    SymbolicTensor::validate_observed_interface(
                        expression,
                        &external,
                        LeafInference::for_declaration(&external),
                    )?;
                    Ok(())
                })
        });
        let head = rebuilt_head
            .filter(|head| !InterfaceInference::intrinsic_normalization_head(*head))
            .or_else(|| {
                source
                    .into_iter()
                    .flat_map(|source| source.get_all_symbols(true))
                    .find(|head| !InterfaceInference::intrinsic_normalization_head(*head))
            })
            .or(rebuilt_head);
        assert!(
            matches!(&checked, Some(Ok(()))),
            "normalization callback contract violated in {operation}; head {head:?}; source {source:?}; result {expression}; admitted interface {structure:?}: {checked:?}"
        );
    }
}

impl<'a> BorrowedTensorParts<'a> {
    fn is_scalar(&self) -> bool {
        self.structure.canonical().is_scalar()
    }

    /// Borrow published construction history. Raw mutation requires an explicit
    /// publication through `finish_parts`; algebra never re-admits its inputs.
    pub(super) fn explicit_index_occurrences(
        &self,
    ) -> InferenceResult<&'a ExplicitIndexOccurrences> {
        self.proofs
            .occurrences
            .get()
            .map(Arc::as_ref)
            .ok_or_else(|| {
                TensorInferenceError::invalid("tensor parts must be published before algebra; finish mutated parts explicitly")
            })
    }

    /// Attach a result of established algebra without making the normalized
    /// payload a second authority for its tensor metadata.
    fn with_established_expression(
        &self,
        expression: Atom,
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        self.explicit_index_occurrences()?;
        let mut proofs = self.proofs.clone();
        proofs.publish_interfaces(expression.as_view(), self.structure.clone(), []);
        Ok(SymbolicTensor {
            is_metric: matches!(expression.as_view(), AtomView::Fun(fun) if fun.get_symbol() == ETS.metric),
            is_composite: !matches!(expression.as_view(), AtomView::Fun(_)),
            expression,
            structure: self.structure.clone(),
            work: Default::default(),
            proofs,
        })
    }

    pub(super) fn normalization_is_intrinsic(&self) -> bool {
        let check = || {
            !self.expression.as_view().needs_normalization()
                && InterfaceInference::normalization_is_intrinsic(self.expression.as_view())
        };
        self.work
            .map_or_else(check, |work| *work.intrinsic.get_or_init(check))
    }

    pub(super) fn product_is_normalized(&self) -> bool {
        self.proofs.is_admitted()
            && !self.structure.has_unresolved_ports()
            && *self.proofs.normalized_product.get_or_init(|| {
                !self.expression.as_view().needs_normalization()
                    && !SymbolicTensor::<PartialStructure>::has_unresolved_ports(
                        self.expression.as_view(),
                    )
                    && matches!(
                        SymbolicTensor::<PartialStructure>::normalize_products(
                            self.expression.as_view()
                        ),
                        AtomOrView::View(_)
                    )
            })
    }

    pub(super) fn established_interface_is_valid(&self) -> bool {
        self.explicit_index_occurrences().is_ok()
    }

    #[cfg(test)]
    pub(super) fn algebra_preserves_interface(&self) -> bool {
        self.structure.open_positions().is_empty()
            && !self.expression.as_view().needs_normalization()
            && InterfaceInference::default()
                .algebra_preserves_leaf_interfaces(self.expression.as_view())
            && self.established_interface_is_valid()
    }

    pub(super) fn rewrites_preserve_interface(&self) -> bool {
        let check = || {
            self.structure.open_positions().is_empty()
                && !self.expression.as_view().needs_normalization()
                && InterfaceInference::default()
                    .rewrites_preserve_leaf_interfaces(self.expression.as_view())
                && self.established_interface_is_valid()
        };
        self.work
            .map_or_else(check, |work| *work.rewrite.get_or_init(check))
    }

    pub(super) fn with_checked_expression(
        &self,
        expression: Atom,
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        if expression.is_zero() && expression.as_view() != self.expression.as_view() {
            return self.with_established_expression(expression);
        }
        if expression.as_view() != self.expression.as_view()
            && !self.structure.open_positions().is_empty()
        {
            return self.finish_checked_open_expression(expression);
        }
        self.finish_checked_expression(expression, self.structure, LeafInference::Observe)
    }

    fn finish_checked_open_expression(
        &self,
        expression: Atom,
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        // The source already supplies the interface; recover only its encoded
        // boundary identities here. Validate scopes on the actual result once.
        let encoded = InterfaceInference {
            leaf_inference: LeafInference::ObserveEncoded,
            ..Default::default()
        }
        .infer_view(self.expression.as_view())
        .and_then(|interface| {
            InterfaceInference::merge_explicit_interface_sequence(&[interface])
                .map(|interface| interface.canonicalize_open_ports())
        });
        if let Ok(encoded) = encoded
            && (encoded.open_positions().is_empty()
                || InterfaceInference::unambiguous_anonymous_ports(&encoded))
        {
            // Retain occurrence/axis identities when explicit ports move past
            // AUTO ports. Observing only anonymous positions can silently swap
            // tensor components after a callback changes normalized factor order.
            let partial = PartialStructure::from_logical_slots(
                encoded.logical_slots().into_iter().map(|slot| {
                    slot.rep().slot(match slot.aind {
                        PartialIndex::Explicit(index) => {
                            InterfaceInference::partial_index(index, LeafInference::Observe)
                        }
                        open => open,
                    })
                }),
            );
            if partial.canonical() == self.structure.canonical()
                && (encoded.open_positions().is_empty()
                    || partial.logical_slots() == self.structure.logical_slots())
            {
                return self.finish_checked_expression(
                    expression,
                    &encoded,
                    LeafInference::ObserveEncoded,
                );
            }
        }
        // Closed scalar factors do not identify anonymous tensor ports. Keep
        // the ordered non-scalar factors as witnesses. One changed rank-one
        // anonymous factor has one possible port correspondence if every other
        // factor stays exactly in place; multiple changed factors do not have
        // that certificate, even when their ranks and representations agree.
        self.explicit_index_occurrences()?;
        let tensor_factors =
            |expression: AtomView<'_>| -> InferenceResult<Vec<(Atom, PartialStructure)>> {
                let factors = match expression {
                    AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                    expression => vec![expression],
                };
                let mut inference = InterfaceInference {
                    leaf_inference: LeafInference::Observe,
                    ..Default::default()
                };
                let mut tensors = Vec::new();
                for factor in factors {
                    // Scalar weights carry no anonymous-axis identity. Match
                    // ordinary product inference here; the checked finisher
                    // still validates the complete result's syntax and scopes.
                    if !InterfaceInference::has_structured_syntax(factor) {
                        continue;
                    }
                    let interface = inference.infer_view(factor)?;
                    if !interface.canonical().is_scalar() {
                        tensors.push((factor.to_owned(), interface));
                    }
                }
                Ok(tensors)
            };
        let retained = tensor_factors(self.expression.as_view())?;
        let rewritten = tensor_factors(expression.as_view())?;
        let mut changed_vectors = 0;
        let retained_ports = !retained.is_empty()
            && retained.len() == rewritten.len()
            && retained
                .iter()
                .zip(&rewritten)
                .all(|((before, left), (after, right))| {
                    if before == after {
                        return true;
                    }
                    changed_vectors += 1;
                    changed_vectors == 1
                        && left.canonical().order() == 1
                        && left.open_positions().len() == 1
                        && left.logical_slots() == right.logical_slots()
                });
        if retained_ports {
            return self.finish_checked_expression(
                expression,
                self.structure,
                LeafInference::Observe,
            );
        }
        // Report invalid result syntax before refusing an otherwise valid
        // rewrite whose unresolved identities have no safe witness. This
        // observation cannot admit the rewrite or replay a normalizer.
        SymbolicTensor::<PartialStructure>::validate_interface(&expression, self.structure)?;
        Err(TensorInferenceError::invalid(
            "cannot certify unresolved tensor port identities after transformation",
        ))
    }

    fn finish_checked_expression(
        &self,
        expression: Atom,
        observed: &PartialStructure,
        policy: LeafInference,
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        if expression.as_view() == self.expression.as_view() {
            return Ok(SymbolicTensor {
                work: self.work.cloned().unwrap_or_default(),
                proofs: self.proofs.clone(),
                expression,
                structure: self.structure.clone(),
                is_metric: self.is_metric,
                is_composite: self.is_composite,
            });
        }
        let ObservedTensor {
            interface: inferred,
            occurrences,
            interfaces,
            leaf_interfaces,
        } = SymbolicTensor::<PartialStructure>::validate_observed_parts(
            &expression,
            observed,
            policy,
            ObservationPurpose::Publish,
        )?;
        if policy == LeafInference::ObserveEncoded && !expression.as_view().is_zero() {
            // Explicit ports can move around unresolved ones, but relabelling
            // same-representation AUTO axes by their new syntax positions would
            // change the tensor's component map.
            if !InterfaceInference::encoded_open_orders_match(observed, &inferred) {
                return Err(TensorInferenceError::invalid(
                    "transformed expression changes the ordered identities of unresolved tensor ports",
                ));
            }
        }
        let mut occurrences = match occurrences {
            _ if expression.is_zero() => ExplicitIndexOccurrences::from_interface(self.structure),
            Some(occurrences) => occurrences,
            None => InterfaceInference::scoped_occurrences(
                expression.as_view(),
                &mut SlotMatcher::default(),
            )?,
        };
        occurrences.branch_max(self.explicit_index_occurrences()?);
        occurrences.validate()?;
        let structure = SymbolicTensor::canonical_structure(self.structure.clone());
        let mut proofs = TensorProofs {
            admitted_interfaces: OnceLock::from(Arc::new(interfaces)),
            leaf_interfaces: leaf_interfaces.into(),
            ..TensorProofs::default()
        };
        proofs.publish_interfaces(expression.as_view(), structure.clone(), [self.proofs]);
        SymbolicTensor::from_published_parts(
            expression,
            structure,
            Arc::new(occurrences),
            proofs.admitted_interfaces,
            proofs.leaf_interfaces,
        )
        .finish_parts()
    }

    /// An oriented product may normalize to a power which raw tensor grammar
    /// rejects. Independently verify both copies when that exact power was
    /// published by the source, then check every other arithmetic operand.
    #[cfg(debug_assertions)]
    fn observe_coalesced_products(
        &self,
        expression: &Atom,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        let facts = TensorProofs::default();
        let policy = LeafInference::for_declaration(self.structure);
        let mut pending = vec![expression.as_view()];
        let mut coalesced = false;
        while let Some(value) = pending.pop() {
            match value {
                AtomView::Mul(product) => pending.extend(product.iter()),
                AtomView::Add(sum) => pending.extend(sum.iter()),
                AtomView::Pow(power)
                    if power.get_base_exp().1 == Atom::num(2).as_view()
                        && self
                            .proofs
                            .retained_interface(value)
                            .is_some_and(|interface| {
                                interface
                                    .logical_slots()
                                    .iter()
                                    .any(|port| !port.rep().rep.is_self_dual())
                            }) =>
                {
                    let base = power.get_base_exp().0.to_owned();
                    let observed = SymbolicTensor::observed_interface(
                        &base,
                        policy,
                        ObservationPurpose::Validate,
                    )?
                    .interface;
                    let copies = InterfaceInference::admitted_oriented_square(
                        self.proofs,
                        value,
                        &observed,
                        power.get_base_exp().1,
                    )
                    .ok_or_else(|| {
                        TensorInferenceError::invalid(
                            "oriented power differs from its two admitted copies",
                        )
                    })?;
                    facts.retain_interface(base.as_view(), observed);
                    facts.retain_interface(value, copies);
                    coalesced = true;
                }
                _ => {
                    let observed = SymbolicTensor::observed_interface(
                        &value.to_owned(),
                        policy,
                        ObservationPurpose::Validate,
                    )?
                    .interface;
                    facts.retain_interface(value, observed);
                }
            }
        }
        if !coalesced
            || !self.structure.open_positions().is_empty()
            || !facts
                .admitted_interface(expression.as_view())
                .is_some_and(|interface| {
                    InterfaceInference::additive_interfaces_match(self.structure, &interface)
                })
        {
            return Err(TensorInferenceError::invalid(
                "oriented products do not establish the result interface",
            ));
        }
        InterfaceInference::produced_occurrences(
            expression.as_view(),
            &mut SlotMatcher::default(),
            &facts,
        )
    }

    /// Publish an operation whose owner preserves the admitted tensor meaning.
    /// Debug observation never supplies production metadata.
    pub(crate) fn with_produced_expression(
        &self,
        expression: Atom,
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        self.explicit_index_occurrences()?;
        if expression.as_view() == self.expression.as_view() {
            return Ok(SymbolicTensor {
                work: self.work.cloned().unwrap_or_default(),
                proofs: self.proofs.clone(),
                expression,
                structure: self.structure.clone(),
                is_metric: self.is_metric,
                is_composite: self.is_composite,
            });
        }
        // Callback contracts preserve free indices but may rename internal
        // dummies. Account those identities without re-inferring the interface.
        if !expression.is_zero() && !self.normalization_is_intrinsic() {
            return self.with_materialized_expression("with_produced_expression", expression);
        }
        #[cfg(debug_assertions)]
        {
            #[cfg(test)]
            let _debug = DebugInference::enter();
            // The checked boundary also verifies encoded AUTO identities. Its
            // observed metadata is discarded: the producer's history remains
            // authoritative, including cancelled dummy indices.
            let checked = self
                .with_checked_expression(expression.clone())
                .and_then(|value| value.explicit_index_occurrences().cloned())
                .or_else(|raw_error| {
                    self.observe_coalesced_products(&expression)
                        .or(Err(raw_error))
                });
            assert!(
                checked.is_ok(),
                "tensor producer contract violated: source {}; result {expression}; admitted interface {:?}: {checked:?}",
                self.expression,
                self.structure,
            );
            let checked = checked.expect("verified producer interface");
            let retained = self.explicit_index_occurrences()?;
            let bounded = checked.is_bounded_by(retained)
                || (self.structure.has_unresolved_ports() && {
                    // Typed products omit physical AUTO markers from their public
                    // occurrence history. Allow the exact markers already in the
                    // source bytes, with their original owners and multiplicities.
                    let mut visible = ExplicitIndexOccurrences::from_atom(
                        self.expression.as_view(),
                        &mut SlotMatcher::default(),
                    );
                    visible.retain(|slot| {
                        CookSettings::index_has_open_port(slot.aind()).unwrap_or(false)
                    });
                    visible.branch_max(retained);
                    checked.is_bounded_by(&visible)
                });
            assert!(
                bounded,
                "tensor producer introduced unreported index occurrences: source {}; result {expression}; observed {:?}; retained {:?}",
                self.expression,
                checked.oriented_entries().collect::<Vec<_>>(),
                retained.oriented_entries().collect::<Vec<_>>(),
            );
        }
        self.with_established_expression(expression)
    }

    pub(super) fn with_algebra_result(
        &self,
        expression: Atom,
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        // Scalar algebra can merge or reorder automatic ports, which have no
        // explicit history to compare against, for example `p(mink(4))^2`
        // after expansion. Check such a changed result instead of trusting it.
        if expression.as_view() != self.expression.as_view()
            && !expression.is_zero()
            && !self.structure.open_positions().is_empty()
        {
            return self.with_checked_expression(expression);
        }
        self.with_produced_expression(expression)
    }

    pub(super) fn with_rewritten_expression(
        &self,
        expression: Atom,
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        self.with_produced_expression(expression)
    }

    pub(crate) fn with_materialized_expression(
        &self,
        operation: &str,
        expression: Atom,
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        if expression.as_view() == self.expression.as_view() {
            return self.with_produced_expression(expression);
        }
        self.with_materialized_history(
            operation,
            expression,
            self.structure.clone(),
            self.explicit_index_occurrences()?,
        )
    }

    /// An intrinsic rewrite can publish its generated leaves and promised
    /// boundary together. A callback's changed result must be checked first.
    pub(crate) fn with_generated_expression(
        &self,
        operation: &str,
        expression: Atom,
        generated: &[TensorProofs],
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        let mut proofs = self.proofs.clone();
        let unchanged = expression.as_view() == self.expression.as_view();
        if unchanged && generated.is_empty() {
            return BorrowedTensorParts {
                proofs: &proofs,
                ..*self
            }
            .with_produced_expression(expression);
        }
        if !unchanged && (expression.is_zero() || self.normalization_is_intrinsic()) {
            let structure = SymbolicTensor::canonical_structure(self.structure.clone());
            proofs.publish_interfaces(expression.as_view(), structure.clone(), generated);
            return Self::finish_materialized_parts(
                operation,
                expression,
                structure,
                self.explicit_index_occurrences()?,
                proofs,
            );
        }
        let established = unchanged.then(|| self.structure.clone());
        proofs.publish_interfaces(expression.as_view(), established, generated);
        BorrowedTensorParts {
            proofs: &proofs,
            work: None,
            ..*self
        }
        .with_materialized_expression(operation, expression)
    }

    /// A scalar-argument normalizer may change the identities inside a leaf.
    /// Reuse its admitted argument roles when the same generic head survives.
    fn callback_leaf_ports(&self, expression: &Atom) -> Option<Vec<Atom>> {
        let (AtomView::Fun(source), AtomView::Fun(result)) =
            (self.expression.as_view(), expression.as_view())
        else {
            return None;
        };
        if source.get_symbol() != result.get_symbol()
            || source.get_nargs() != result.get_nargs()
            || InterfaceInference::is_builtin_tensor(source.get_symbol())
        {
            return None;
        }
        let (true, ports) = self.source_proofs.leaf_interfaces.fact(source.as_view())? else {
            return None;
        };
        let mut written = ports.iter().peekable();
        let mut mapped = Vec::with_capacity(ports.len());
        let mut slots = SlotMatcher::default();
        for (before, after) in source.iter().zip(result.iter()) {
            if written.next_if(|port| port.as_view() == before).is_some() {
                if slots.parse::<LibraryRep, AbstractIndex>(after).is_err()
                    && slots.parse_representation::<LibraryRep>(after).is_err()
                {
                    return None;
                }
                mapped.push(after.to_owned());
            } else if before != after
                && !matches!((before, after), (AtomView::Num(_), AtomView::Num(_)))
            {
                return None;
            }
        }
        written.next().is_none().then_some(mapped)
    }

    /// Arithmetic on admitted operands carries their interfaces. A new
    /// function emitted by user code still needs its argument grammar checked.
    fn admitted_callback_result(
        inputs: &TensorProofs,
        generated: &TensorProofs,
        expression: AtomView<'_>,
    ) -> Option<PartialStructure> {
        let mut pending = vec![expression];
        while let Some(value) = pending.pop() {
            if inputs.retained_interface(value).is_some() || inputs.leaf_interfaces.contains(value)
            {
                continue;
            }
            match value {
                AtomView::Num(_) => {}
                AtomView::Var(variable)
                    if !variable.get_symbol().has_tag(&SPENSO_TAG.tensor)
                        && !InterfaceInference::is_chain_placeholder(value) => {}
                AtomView::Mul(product) => pending.extend(product.iter()),
                AtomView::Add(sum) => pending.extend(sum.iter()),
                AtomView::Pow(power) if matches!(power.get_base_exp().1, AtomView::Num(_)) => {
                    pending.push(power.get_base_exp().0);
                }
                AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.chain => {
                    let arguments = function.iter().collect::<Vec<_>>();
                    let mut grammar = InterfaceInference {
                        leaf_inference: LeafInference::ObserveEncoded,
                        ..Default::default()
                    };
                    let interface = grammar
                        .chain_interface_with(&arguments, |_, factor| {
                            Self::admitted_callback_result(inputs, generated, factor).ok_or_else(
                                || {
                                    TensorInferenceError::invalid(
                                        "a callback emitted an unknown chain factor",
                                    )
                                },
                            )
                        })
                        .ok()?;
                    inputs.retain_boundary(value, &interface);
                }
                AtomView::Fun(function) if function.get_symbol().has_tag(&SPENSO_TAG.broadcast) => {
                    if function.get_nargs() != 1 {
                        return None;
                    }
                    let operand = function.iter().next()?;
                    if !InterfaceInference::has_structured_syntax(operand) {
                        return None;
                    }
                    pending.push(operand);
                }
                AtomView::Fun(function)
                    if composition::is_tensor_leaf_head(function.get_symbol())
                        && function.get_symbol() != SPENSO_TAG.dot
                        && !InterfaceInference::is_compound_metric(
                            function,
                            &mut SlotMatcher::default(),
                        ) =>
                {
                    // A leaf producer supplies written argument roles. In
                    // contrast, a composite's parser ports can be synthesized
                    // from the very boundary this validation must check.
                    let (_, ports) = generated
                        .leaf_interfaces
                        .fact(value)
                        .filter(|(tensor, _)| *tensor)?;
                    let mut arguments = function.iter();
                    if !ports
                        .iter()
                        .all(|port| arguments.any(|argument| argument == port.as_view()))
                    {
                        return None;
                    }
                    inputs.retain_leaf_ports(value, ports).ok()?;
                }
                _ => return None,
            }
        }
        inputs.admitted_interface(expression)
    }

    pub(super) fn with_materialized_parts(
        &self,
        operation: &str,
        expression: Atom,
        structure: PartialStructure,
        mut map_slot: impl FnMut(Slot<LibraryRep, AbstractIndex>) -> Slot<LibraryRep, AbstractIndex>,
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        let retained = self
            .explicit_index_occurrences()?
            .mapped(|slot| Ok::<_, TensorInferenceError>(map_slot(slot)))?;
        self.with_materialized_history(operation, expression, structure, &retained)
    }

    fn with_materialized_history(
        &self,
        operation: &str,
        expression: Atom,
        structure: PartialStructure,
        retained: &ExplicitIndexOccurrences,
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        let mut produced = self.proofs.clone();
        if !expression.is_zero()
            && expression.as_view() != self.expression.as_view()
            && !self.normalization_is_intrinsic()
        {
            // User normalizers may erase a free index or replace the entire
            // tensor. Compare the output of the callback before publishing the
            // producer's declared boundary. Existing regional facts cover
            // unchanged operands and arithmetic on them without observation.
            // Generated regional declarations can be this materializer's
            // promise, even when they belong to a draft rooted at the input.
            // Start verification with the independent input admission.
            produced = self.source_proofs.clone();
            produced.publish_interfaces(expression.as_view(), None, []);
            if let Some(ports) = self.callback_leaf_ports(&expression) {
                produced.retain_leaf_ports(expression.as_view(), ports)?;
            }
            let external =
                InterfaceInference::merge_explicit_interface_sequence(slice::from_ref(&structure))?;
            if let Some(actual) =
                Self::admitted_callback_result(&produced, self.proofs, expression.as_view())
            {
                if !InterfaceInference::additive_interfaces_match(&external, &actual) {
                    return Err(TensorInferenceError::invalid(
                        "transformed expression does not preserve a compatible tensor interface",
                    ));
                }
                if std::ptr::eq(self.source_proofs, self.proofs) {
                    let structure = SymbolicTensor::canonical_structure(structure);
                    produced.set_projected_boundary(expression.as_view(), structure.clone());
                    return Self::finish_materialized_parts(
                        operation, expression, structure, retained, produced,
                    );
                }
            } else {
                // The user callback emitted previously unknown syntax. Admit
                // its output once without invoking its normalizer a second time.
                let observed = SymbolicTensor::<PartialStructure>::from_observed_parts(
                    expression.clone(),
                    external.clone(),
                    LeafInference::for_declaration(&external),
                )?;
                produced = observed.proofs;
            }
        }
        let structure = SymbolicTensor::canonical_structure(structure);
        produced.publish_interfaces(expression.as_view(), structure.clone(), [self.proofs]);
        Self::finish_materialized_parts(operation, expression, structure, retained, produced)
    }

    fn finish_materialized_parts(
        operation: &str,
        expression: Atom,
        structure: PartialStructure,
        retained: &ExplicitIndexOccurrences,
        produced: TensorProofs,
    ) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        let mut occurrences = InterfaceInference::produced_occurrences(
            expression.as_view(),
            &mut SlotMatcher::default(),
            &produced,
        )?;
        occurrences.branch_max(retained);
        occurrences.branch_max(&ExplicitIndexOccurrences::from_interface(&structure));
        occurrences.validate()?;
        #[cfg(debug_assertions)]
        SymbolicTensor::<PartialStructure>::assert_produced_interface(
            operation,
            &expression,
            &structure,
        );
        #[cfg(not(debug_assertions))]
        let _ = operation;
        Ok(SymbolicTensor::from_published_parts(
            expression,
            structure,
            Arc::new(occurrences),
            produced.admitted_interfaces,
            produced.leaf_interfaces,
        ))
    }

    pub(super) fn with_transformed_expression(
        &self,
        expression: Atom,
    ) -> InferenceResult<(SymbolicTensor<PartialStructure>, bool)> {
        if expression.as_view() == self.expression.as_view() {
            return Ok((
                SymbolicTensor {
                    work: self.work.cloned().unwrap_or_default(),
                    proofs: self.proofs.clone(),
                    expression,
                    structure: self.structure.clone(),
                    is_metric: self.is_metric,
                    is_composite: self.is_composite,
                },
                false,
            ));
        }
        if expression.is_zero() {
            return Ok((self.with_established_expression(expression)?, false));
        }
        if self.is_scalar() || InterfaceInference::has_structured_syntax(expression.as_view()) {
            return SymbolicTensor::<PartialStructure>::infer_with_contractions(expression);
        }
        Err(TensorInferenceError::invalid(
            "Tensor transformation removed the tensor syntax of a non-scalar expression",
        ))
    }
}

impl SymbolicTensor<PartialStructure> {
    pub(super) fn borrowed_parts(&self) -> BorrowedTensorParts<'_> {
        BorrowedTensorParts {
            work: Some(&self.work),
            ..self.borrowed_with_proofs(&self.proofs)
        }
    }

    /// Borrow the established payload with an operation's independent fact cache.
    /// The supplied proofs must retain this tensor's interface and index history.
    pub(crate) fn borrowed_with_proofs<'a>(
        &'a self,
        proofs: &'a TensorProofs,
    ) -> BorrowedTensorParts<'a> {
        BorrowedTensorParts {
            expression: &self.expression,
            structure: &self.structure,
            is_metric: self.is_metric,
            is_composite: self.is_composite,
            work: None,
            proofs,
            source_proofs: &self.proofs,
        }
    }

    pub(crate) fn explicit_index_occurrences(&self) -> InferenceResult<&ExplicitIndexOccurrences> {
        self.borrowed_parts().explicit_index_occurrences()
    }

    /// Iterate retained explicit indices, including identities hidden by cancellation.
    /// Transporters use these together with the logical interface; see
    /// [`Self::with_renamed_indices`] for a metadata-preserving relabeling.
    ///
    /// # Errors
    /// Returns an admission error if the tensor's retained metadata is invalid.
    pub fn explicit_indices(&self) -> InferenceResult<impl Iterator<Item = AbstractIndex> + '_> {
        Ok(self
            .explicit_index_occurrences()?
            .entries()
            .map(|(slot, _)| slot.aind()))
    }

    /// Construct parts whose producer has established their complete tensor interface.
    /// This performs no inference, validation, or normalization. Callers must
    /// establish compatible ports in every summand, valid tensor signatures,
    /// placeholder scopes and index multiplicities, and account for callback
    /// effects. Preserve the chosen logical order, including the rank of zero.
    /// Raw external expressions belong at `from_declared_parts` instead.
    /// The producer supplies conservative construction history, including dummy
    /// uses hidden by normalization. Publication records exact payload keys but
    /// never discovers tensor ports from the payload.
    pub(crate) fn from_validated_parts(
        expression: Atom,
        structure: PartialStructure,
        occurrences: Arc<ExplicitIndexOccurrences>,
    ) -> Self {
        let structure = Self::canonical_structure(structure);
        let interfaces = AdmittedInterfaces::from_region(expression.as_view(), structure.clone());
        Self::from_published_parts(
            expression,
            structure,
            occurrences,
            OnceLock::from(Arc::new(interfaces)),
            Default::default(),
        )
    }

    pub(super) fn canonical_structure(structure: PartialStructure) -> PartialStructure {
        if structure.has_unresolved_ports() {
            structure.canonicalize_open_ports()
        } else {
            structure
        }
    }

    /// The supplied maps already belong to this exact expression and carry its
    /// canonical root. Only the new index history and validation are certified.
    pub(super) fn from_published_parts(
        expression: Atom,
        structure: PartialStructure,
        occurrences: Arc<ExplicitIndexOccurrences>,
        admitted_interfaces: OnceLock<Arc<AdmittedInterfaces>>,
        leaf_interfaces: RawLeafInterfaces,
    ) -> Self {
        let is_metric =
            matches!(expression.as_view(), AtomView::Fun(fun) if fun.get_symbol() == ETS.metric);
        let is_composite = !matches!(expression.as_view(), AtomView::Fun(_));
        Self {
            expression,
            structure,
            is_metric,
            is_composite,
            work: Default::default(),
            proofs: TensorProofs {
                occurrences: OnceLock::from(occurrences),
                admitted_interfaces,
                leaf_interfaces,
                ..TensorProofs::default()
            },
        }
    }

    pub(crate) fn normalization_is_intrinsic(&self) -> bool {
        self.borrowed_parts().normalization_is_intrinsic()
    }

    /// Publish an identity implemented by our algebra, preserving its proven
    /// logical interface. Normalization callbacks obey the same free-index
    /// preservation contract; debug checks never supply production metadata.
    pub(crate) fn with_identity_result(
        &self,
        expression: Atom,
        observations: Option<DomainObservations>,
    ) -> InferenceResult<Self> {
        self.explicit_index_occurrences()?;
        let unchanged = expression == self.expression;
        if unchanged && observations.is_none() {
            return Ok(self.clone());
        }
        // Free-index preservation does not certify identity freedom or dot
        // notation: a callback can introduce closed tensor work. New producer
        // facts describe a changed payload only when callbacks were excluded.
        let preserve_observations = unchanged || self.normalization_is_intrinsic();
        let mut result = if unchanged {
            // Even a syntactic fixed point can acquire stronger producer facts.
            // Preserve its construction history and all payload-specific proofs.
            self.clone()
        } else if expression.is_zero()
            || (preserve_observations
                && observations
                    .as_ref()
                    .is_some_and(DomainObservations::is_terminal_polynomial))
        {
            // Terminal certificates exclude internal connections. Every visible
            // index is already part of the retained interface/history.
            self.with_produced_expression(expression)?
        } else {
            // General color/gamma identities can introduce internal contractions.
            // Account for those occurrences using the producer's regional facts.
            self.with_materialized_expression("with_identity_result", expression)?
        };
        if preserve_observations && let Some(observations) = observations {
            // Keep equivalent or stronger facts, including their region and
            // scope caches, when an unchanged producer republishes its proof.
            if !result.work.observations().is_some_and(|previous| {
                previous.is_terminal_polynomial()
                    && observations.is_terminal_polynomial()
                    && (!previous.needs_dot_normalization()
                        || observations.needs_dot_normalization())
            }) {
                let observations = if unchanged && let Some(previous) = result.work.observations() {
                    observations.retain_regions(previous)
                } else {
                    observations
                };
                result.work.set_observations(observations);
            }
        }
        if self.work.observations().is_some() {
            result.observe_replacement(self);
        }
        Ok(result)
    }

    /// Publish a built-in materializer's output, retaining this logical interface
    /// and every previously reserved index. Newly emitted indices are accounted
    /// for without inferring the tensor interface. See the
    /// [producer transport example](Self::with_produced_parts).
    ///
    /// # Errors
    /// Rejects unpublished source metadata or invalid index multiplicities.
    ///
    /// # Panics
    /// Panics in debug builds if the materializer changes the admitted interface.
    pub fn with_materialized_expression(
        &self,
        operation: &str,
        expression: Atom,
    ) -> InferenceResult<Self> {
        self.borrowed_parts()
            .with_materialized_expression(operation, expression)
    }

    /// Finish a producer's private, already scoped output. Its index history
    /// still belongs to the source until this materializer accounts for new uses.
    pub(crate) fn with_published_materialization(
        &self,
        operation: &str,
        expression: Atom,
        mut proofs: TensorProofs,
    ) -> InferenceResult<Self> {
        if expression.as_view() == self.expression.as_view() {
            if !expression.is_zero()
                && !proofs.owns_published_interface(expression.as_view(), &self.structure)
            {
                proofs.publish_interfaces(expression.as_view(), self.structure.clone(), []);
            }
            return self
                .borrowed_with_proofs(&proofs)
                .with_produced_expression(expression);
        }
        let source = self.borrowed_with_proofs(&proofs);
        let structure = Self::canonical_structure(self.structure.clone());
        if !expression.is_zero()
            && (!source.normalization_is_intrinsic()
                || !proofs.owns_published_interface(expression.as_view(), &structure))
        {
            return source.with_materialized_expression(operation, expression);
        }
        let retained = Arc::clone(proofs.occurrences.get().ok_or_else(|| {
            TensorInferenceError::invalid("published materializer is missing construction history")
        })?);
        BorrowedTensorParts::finish_materialized_parts(
            operation, expression, structure, &retained, proofs,
        )
    }

    /// Publish a materializer which also transports the logical interface.
    /// The slot map transforms all retained history, including erased dummies;
    /// emitted index occurrences add newly allocated identities. Only this index
    /// accounting reads the output in release builds, never interface inference.
    /// See the [producer transport example](Self::with_produced_parts).
    ///
    /// # Errors
    /// Rejects unpublished metadata or invalid index multiplicities.
    ///
    /// # Panics
    /// Panics in debug builds if the payload disagrees with its supplied interface.
    pub fn with_materialized_parts(
        &self,
        operation: &str,
        expression: Atom,
        structure: PartialStructure,
        map_slot: impl FnMut(Slot<LibraryRep, AbstractIndex>) -> Slot<LibraryRep, AbstractIndex>,
    ) -> InferenceResult<Self> {
        self.borrowed_parts()
            .with_materialized_parts(operation, expression, structure, map_slot)
    }

    /// Inspect the exact raw leaf facts retained by this carrier.
    #[cfg(test)]
    pub(crate) fn admitted_leaf_interfaces(&self) -> LeafInterfaceCache {
        Arc::clone(self.proofs.leaf_interfaces.get_or_init(Default::default))
    }

    /// Cached facts belong to this exact payload/interface pair. Unchecked
    /// storage construction starts with no facts; mutable trait access clears them.
    pub(super) fn established_interface_is_valid(&self) -> bool {
        self.borrowed_parts().established_interface_is_valid()
    }

    pub(super) fn rewrites_preserve_interface(&self) -> bool {
        self.borrowed_parts().rewrites_preserve_interface()
    }

    /// Construct an unresolved tensor from a named, canonically ordered signature.
    /// The signature supplies both logical axes and written parser ports.
    /// See the [materialization example](Self::materialized).
    ///
    /// # Errors
    /// Rejects an unnamed signature or a user normalizer that changes its interface.
    pub fn from_signature(
        value: &Canonicalized<ExplicitKey<AbstractIndex>>,
    ) -> InferenceResult<Self> {
        let signature = value.canonical();
        let name = signature
            .global_name
            .ok_or_else(|| TensorInferenceError::invalid("tensor structure has no name"))?;
        let arguments = signature.additional_args.as_deref().unwrap_or_default();
        let canonical = signature.external_reps_iter().collect::<Vec<_>>();
        let logical = value.layout().canonical_to_logical(&canonical);
        let axes = (0..logical.len()).collect::<Vec<_>>();
        let canonical_axes = value.layout().logical_to_canonical(&axes);
        let owner = AbstractIndex::fresh_open_owner();
        let ports = canonical
            .iter()
            .zip(canonical_axes)
            .map(|(representation, axis)| {
                representation.slot::<AbstractIndex, _>(AbstractIndex::Open { owner, axis })
            })
            .collect::<Vec<_>>();
        let raw_ports = ports
            .iter()
            .map(IsAbstractSlot::to_atom)
            .collect::<Vec<_>>();
        let expression = FunctionBuilder::new(name)
            .add_args(arguments)
            .add_args(&raw_ports)
            .finish();
        let storage = OrderedStructure::new(ports).into_canonical();
        let proofs = SymbolicTensor::<OrderedStructure>::constructor_proofs(
            &expression,
            &storage,
            name,
            raw_ports,
        )?;
        let structure =
            PartialStructure::from_logical_slots(logical.into_iter().enumerate().map(
                |(position, representation)| representation.slot(PartialIndex::open(position)),
            ));
        let mut value = Self::from_validated_parts(
            expression,
            structure,
            Arc::clone(
                proofs
                    .occurrences
                    .get()
                    .expect("the signature publishes construction history"),
            ),
        );
        value.proofs.inherit_interfaces([&proofs]);
        let value = value.finish_parts()?;
        #[cfg(debug_assertions)]
        Self::assert_constructor_interface("from_signature", &value.expression, &value.structure);
        Ok(value)
    }

    /// Infer an ordered interface from symbolic syntax, lowering unresolved powers.
    pub fn infer(atom: Atom) -> InferenceResult<Self> {
        #[cfg(feature = "reference-cases")]
        let _phase =
            crate::reference_cases::timing::scope(crate::reference_cases::timing::Phase::Admission);
        Self::infer_with_contractions(atom).map(|(value, _)| value)
    }

    fn infer_with_contractions(atom: Atom) -> InferenceResult<(Self, bool)> {
        // Tensor-free input has no tensor interface or index-use history.
        // This does not certify normalization or any downstream algebra guarantee.
        if !InterfaceInference::contains_tensor_symbol(atom.as_view()) {
            let scalar = PartialStructure::from_logical_slots([]);
            let value = Self::from_validated_parts(
                atom,
                scalar.clone(),
                Arc::new(ExplicitIndexOccurrences::default()),
            );
            // Arithmetic can expose an opaque coefficient of this tensor-free
            // input as a separate graph leaf, for example when a sign flattens
            // a product. Every such coefficient already passed scalar admission.
            value.expression.visitor(&mut |node| {
                if matches!(node, AtomView::Fun(_)) {
                    value.proofs.retain_interface(node, scalar.clone());
                    return false;
                }
                true
            });
            return Ok((value, false));
        }
        Self::admit_with_contractions(atom)
    }

    fn admit_with_contractions(atom: Atom) -> InferenceResult<(Self, bool)> {
        let atom = DotNormalizer::nested_vectors(atom.as_view());
        let atom = InterfaceInference::lower_tensor_powers(atom.as_view())?.unwrap_or(atom);
        let mut inference = InterfaceInference {
            leaf_inference: LeafInference::Observe,
            purpose: ObservationPurpose::Publish,
            ..InterfaceInference::default()
        };
        let (structure, occurrences) = inference.infer_observed(atom.as_view())?;
        let raw_rank = structure.canonical().order();
        let structure = InterfaceInference::merge_explicit_interface_sequence(&[structure])?
            .canonicalize_open_ports();
        let contracted = structure.canonical().order() != raw_rank;
        let (atom, unchanged) = match Self::normalize_products(atom.as_view()) {
            AtomOrView::View(_) => (atom, true),
            normalized => (normalized.into_owned(), false),
        };
        let value = inference.publish_admission(atom, structure, occurrences);
        // Admission publishes its observed interface and conservative history
        // before normalization. Payload-specific optimization facts additionally
        // require unchanged, intrinsic syntax.
        let intrinsic = !value.expression.as_view().needs_normalization()
            && inference.normalization_is_intrinsic_cached(value.expression.as_view());
        let _ = value.work.intrinsic.set(intrinsic);
        if unchanged && intrinsic {
            let _ = value
                .proofs
                .normalized_product
                .set(!Self::has_unresolved_ports(value.expression.as_view()));
            if value.structure.open_positions().is_empty()
                && inference
                    .algebra_preserves_leaf_interfaces_impl(value.expression.as_view(), true)
            {
                let _ = value.work.rewrite.set(true);
            }
        }
        Ok((value, contracted))
    }

    /// Raise a typed value to a scalar power without losing its logical interface.
    /// Non-scalar powers use the existing tensor-product plan for each occurrence.
    /// See the [example](SymbolicTensor).
    ///
    /// # Errors
    ///
    /// Rejects invalid operands, a tensor-valued exponent, unsupported non-scalar powers,
    /// or an index occurring more than twice in the resulting multiplicative scope.
    pub fn pow(&self, exponent: &Self) -> InferenceResult<Self> {
        if !self.established_interface_is_valid() || !exponent.established_interface_is_valid() {
            return Err(TensorInferenceError::invalid(
                "power operands have invalid tensor interfaces",
            ));
        }
        if !exponent.is_scalar() {
            return Err(TensorInferenceError::invalid(
                "a tensor exponent must be scalar",
            ));
        }
        let exponent_atom = Self::scalar_power_expression(
            exponent.expression.clone(),
            exponent.has_indexed_power_base(),
        );
        if exponent_atom.is_one() {
            return Ok(self.clone());
        }
        if exponent_atom.is_zero() {
            let result = Self::from_validated_parts(
                Atom::one(),
                PartialStructure::from_logical_slots([]),
                Arc::new(ExplicitIndexOccurrences::default()),
            );
            let _ = result.proofs.normalized_product.set(true);
            return Ok(result);
        }
        let interface =
            InterfaceInference::power_interface(&self.structure, exponent_atom.as_view(), || {
                format!("{}^({})", self.expression, exponent_atom)
            })?;
        if self.is_scalar() {
            let expression = Self::scalar_power_expression(
                self.expression.clone(),
                self.has_indexed_power_base(),
            )
            .pow(&exponent_atom);
            // A scalar base and its exponent own independent closed scopes.
            // Their internal dummies cannot meet factors outside this power.
            self.explicit_index_occurrences()?;
            exponent.explicit_index_occurrences()?;
            let mut result = Self::from_validated_parts(
                expression,
                interface,
                Arc::new(ExplicitIndexOccurrences::default()),
            );
            result
                .proofs
                .inherit_interfaces([&self.proofs, &exponent.proofs]);
            return Ok(result);
        }
        let rational = Rational::try_from(exponent_atom.as_view()).map_err(|_| {
            TensorInferenceError::invalid("a non-scalar tensor power requires an integer exponent")
        })?;
        let negative = rational.numerator().is_negative();
        if rational.denominator() != 1 || (negative && rational.numerator() % 2 != 0) {
            return Err(TensorInferenceError::invalid(
                "a non-scalar tensor power requires a non-negative or even integer exponent",
            ));
        }
        let repetitions = usize::try_from(rational.numerator().abs()).map_err(|_| {
            TensorInferenceError::invalid("tensor power exponent is too large to materialize")
        })?;
        if negative {
            // A negative even power is the reciprocal of the contracted positive
            // power, which `power_interface` has already found to be scalar.
            // Multiplication contracts every pair of copies, explicit ones
            // included, so the inverse writes no label the surrounding scope
            // could meet, exactly as a typed square raised to -1.
            let base = self.with_automatic_ports()?;
            let mut positive = base.clone();
            for _ in 1..repetitions {
                positive = positive.multiply(&base)?;
            }
            let minus_one =
                Self::from_produced_parts(Atom::num(-1), PartialStructure::from_logical_slots([]))?;
            return positive.pow(&minus_one);
        }
        if repetitions == 0 {
            return Self::from_produced_parts(Atom::one(), interface);
        }
        if repetitions == 1 {
            return Ok(self.clone());
        }
        if self.structure.open_positions().is_empty() {
            // Explicit tensor powers are repeated tensor products. Establish
            // their index uses from the admitted base, before Atom coalescence.
            let mut occurrences = ExplicitIndexOccurrences::default();
            let source = self.explicit_index_occurrences()?;
            for _ in 0..repetitions {
                occurrences.append(source);
            }
            occurrences.validate()?;
            let expression = self.expression.pow(&exponent_atom);
            let mut result =
                Self::from_validated_parts(expression, interface, Arc::new(occurrences));
            result
                .proofs
                .inherit_interfaces([&self.proofs, &exponent.proofs]);
            let _ = result
                .proofs
                .normalized_product
                .set(self.borrowed_parts().product_is_normalized());
            return Ok(result);
        }
        let mut result = self.clone();
        for _ in 1..repetitions {
            result = result.multiply(self)?;
        }
        Ok(result)
    }

    /// The base of a self-contracted inverse power with its explicit ports
    /// written as automatic ones. Those labels are internal to the inverse, so
    /// multiplication may give each pair of copies a fresh dummy, and no label
    /// of the inverse meets the surrounding scope or its canonicalization.
    /// Only the open ports of a rank-one base are rewritten: closed scalars
    /// inside it keep their pairs, only tensor and composite heads carry
    /// indices, and other calls keep their index-shaped data.
    fn with_automatic_ports(&self) -> InferenceResult<Self> {
        let explicit = self
            .structure
            .logical_slots()
            .into_iter()
            .filter(|slot| matches!(slot.aind, PartialIndex::Explicit(_)))
            .collect::<Vec<_>>();
        // Named axes of a higher-rank base pair by label, also across terms that
        // order them differently, so its multiplication keeps the labels.
        if explicit.is_empty() || self.rank() != 1 {
            return Ok(self.clone());
        }
        let mut slots = SlotMatcher::default();
        let expression = self.expression.replace_map(|value, _, output| {
            if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(value) {
                **output = if explicit.iter().any(|port| {
                    port.aind == PartialIndex::Explicit(slot.aind()) && port.rep() == slot.rep()
                }) {
                    slot.rep().to_symbolic([])
                } else {
                    value.to_owned()
                };
            } else if Self::has_closed_tensor_scope(value) {
                // A closed scalar keeps its own pairs, even under the port's label.
                **output = value.to_owned();
            } else if let AtomView::Fun(function) = value {
                let symbol = function.get_symbol();
                if symbol.is_scalar()
                    || !(symbol.has_tag(&SPENSO_TAG.tensor)
                        || composition::is_composite_head(symbol, &SPENSO_TAG))
                {
                    **output = value.to_owned();
                }
            }
        });
        Self::infer(expression)
    }

    /// Whether `value` is tensor syntax with no logical port, explicit ones
    /// included: the canonical structure alone omits explicit labels.
    fn has_closed_tensor_scope(value: AtomView<'_>) -> bool {
        InterfaceInference::has_structured_syntax(value)
            && InterfaceInference::replacement_interface(value)
                .is_ok_and(|interface| interface.logical_slots().is_empty())
    }

    /// Divide by a scalar tensor, retaining any scoped contractions in its inverse.
    pub fn divide(&self, denominator: &Self) -> InferenceResult<Self> {
        if !denominator.is_scalar() {
            return Err(TensorInferenceError::invalid(
                "cannot divide by a non-scalar tensor; tensor division requires a scalar denominator",
            ));
        }
        let exponent = Self::from_validated_parts(
            Atom::num(-1),
            PartialStructure::from_logical_slots([]),
            Arc::new(ExplicitIndexOccurrences::default()),
        );
        let inverse = denominator.pow(&exponent)?;
        Ok(self.multiply(&inverse)?)
    }

    /// Finish parts whose producer guarantees payload/interface agreement.
    ///
    /// The caller must establish that the payload agrees with its logical ports.
    /// External expressions and metadata must use [`Self::from_declared_parts`].
    /// Debug builds check payload/interface agreement; release builds trust the
    /// producer and only account for newly published index occurrences.
    /// Zero retains the declared rank. See also the [example](SymbolicTensor).
    ///
    /// # Errors
    ///
    /// Rejects invalid placeholder scopes or explicit-index multiplicities.
    ///
    /// # Panics
    ///
    /// In debug builds, panics if the finished payload disagrees with its interface.
    ///
    /// # Example
    ///
    /// ```
    /// use idenso::tensor::SymbolicTensor;
    /// use spenso::structure::partial::{PartialStructure, PartialStructureExt};
    /// use symbolica::atom::{Atom, AtomCore};
    /// let scalar = SymbolicTensor::from_produced_parts(
    ///     Atom::num(3), PartialStructure::from_logical_slots([]),
    /// )?;
    /// assert!(scalar.is_scalar());
    /// # Ok::<(), idenso::tensor::inference::TensorInferenceError>(())
    /// ```
    pub fn from_produced_parts(
        expression: Atom,
        structure: PartialStructure,
    ) -> InferenceResult<Self> {
        // Constructor-owned ports are authoritative. Nested scopes must already
        // have facts or be composed from known tensor constructors.
        let mut value = Self::from_validated_parts(
            expression,
            structure,
            Arc::new(ExplicitIndexOccurrences::default()),
        );
        let mut occurrences = InterfaceInference::constructor_occurrences(
            value.expression.as_view(),
            &mut SlotMatcher::default(),
            &value.proofs,
            true,
        )?;
        occurrences.branch_max(&ExplicitIndexOccurrences::from_interface(&value.structure));
        occurrences.validate()?;
        value.proofs.occurrences = OnceLock::from(Arc::new(occurrences));
        let value = value.finish_parts()?;
        #[cfg(debug_assertions)]
        Self::assert_constructor_interface(
            "from_produced_parts",
            &value.expression,
            &value.structure,
        );
        Ok(value)
    }

    /// Admit an external expression with its declared ordered interface.
    ///
    /// Unlike attaching parts whose interface was established by an algebra
    /// owner, this checks the payload, tensor scopes, index multiplicities and
    /// declared ports together. Observation does not replay user normalizers.
    /// Zero retains the declared rank and logical order.
    /// See the [example](SymbolicTensor).
    ///
    /// # Errors
    ///
    /// Rejects invalid expression scopes, invalid index multiplicities, or a declared
    /// logical interface that disagrees with the expression.
    pub fn from_declared_parts(
        expression: Atom,
        structure: PartialStructure,
    ) -> InferenceResult<Self> {
        Self::from_observed_parts(expression, structure, LeafInference::Observe)
    }

    /// Observe an unknown user transformation without replaying normalization.
    /// Encoded observation keeps canonicalizer-owned AUTO identities literal.
    pub(super) fn from_observed_parts(
        expression: Atom,
        structure: PartialStructure,
        policy: LeafInference,
    ) -> InferenceResult<Self> {
        let observed = Self::validate_observed_parts(
            &expression,
            &structure,
            policy,
            ObservationPurpose::Publish,
        )?;
        Self::from_observation(expression, structure, observed)
    }

    /// Admit a raw algebra adapter's literal syntax without rewriting its powers or notation.
    pub(crate) fn observe_raw_expression(expression: Atom) -> InferenceResult<Self> {
        let observed = Self::observed_interface(
            &expression,
            LeafInference::ObserveEncoded,
            ObservationPurpose::Publish,
        )?;
        Self::from_observation(expression, observed.interface.clone(), observed)
    }

    fn from_observation(
        expression: Atom,
        structure: PartialStructure,
        observed: ObservedTensor,
    ) -> InferenceResult<Self> {
        let ObservedTensor {
            occurrences,
            mut interfaces,
            leaf_interfaces,
            ..
        } = observed;
        let mut occurrences = match occurrences {
            Some(occurrences) => occurrences,
            None => InterfaceInference::scoped_occurrences(
                expression.as_view(),
                &mut SlotMatcher::default(),
            )?,
        };
        occurrences.branch_max(&ExplicitIndexOccurrences::from_interface(&structure));
        occurrences.validate()?;
        let structure = Self::canonical_structure(structure);
        let leaf_interfaces =
            interfaces.declare_root(expression.as_view(), structure.clone(), leaf_interfaces);
        let value = Self::from_published_parts(
            expression,
            structure,
            Arc::new(occurrences),
            OnceLock::from(Arc::new(interfaces)),
            leaf_interfaces,
        );
        value.borrowed_parts().product_is_normalized();
        Ok(value)
    }

    /// Finish an existing carrier without detaching its established facts.
    /// Unknown parts still undergo scope and multiplicity validation; changed
    /// product syntax or logical layout cannot inherit the original proofs.
    /// See the [example](SymbolicTensor).
    ///
    /// # Errors
    ///
    /// Rejects invalid source metadata, unmatched interface ports, or a contraction
    /// that violates explicit-index multiplicity.
    pub fn finish_parts(mut self) -> InferenceResult<Self> {
        if !self.proofs.is_admitted() {
            return Self::from_declared_parts(self.expression, self.structure)?.finish_parts();
        }
        self.explicit_index_occurrences()?;
        let normalized = if self.proofs.normalized_product.get() == Some(&true) {
            None
        } else {
            match Self::normalize_products(self.expression.as_view()) {
                AtomOrView::View(_) => None,
                normalized => Some(normalized.into_owned()),
            }
        };
        if let Some(expression) = normalized {
            let mut proofs = self.proofs.clone();
            proofs.normalized_product = Default::default();
            self.expression = expression;
            self.is_metric = matches!(self.expression.as_view(), AtomView::Fun(fun) if fun.get_symbol() == ETS.metric);
            self.is_composite = !matches!(self.expression.as_view(), AtomView::Fun(_));
            proofs.publish_interfaces(self.expression.as_view(), self.structure.clone(), []);
            self.replace_admission(proofs);
        }
        // The construction owner already established this interface. Finishing
        // payload spelling cannot contract two declared ports inside one factor,
        // including the nominal ports of zero. Parser inference and explicit
        // composition own their respective contractions and chain closure.
        Ok(self)
    }

    /// Check an arbitrary result against this interface and finish it once.
    /// Observation checks tensor scopes before index multiplicity; scalar results
    /// still check those scopes after their interface has been accepted.
    pub fn with_checked_expression(&self, expression: Atom) -> InferenceResult<Self> {
        self.borrowed_parts().with_checked_expression(expression)
    }

    /// Build a symbolic evaluator without expanding or contracting this payload.
    /// Component-valued tensors use the network/library materialization boundary;
    /// symbolic tensor leaves may instead be explicit evaluator inputs.
    pub fn evaluator<P: AtomCore>(&self, parameters: &[P]) -> EvaluatorBuilder<'_> {
        self.expression.evaluator(parameters)
    }

    /// A plain scalar symbol absent from `payload` and `selected`, which marks
    /// selected occurrences during one expansion. Each candidate is requested
    /// without attributes through the fallible builder: a name which another
    /// owner registered differently is refused and skipped instead of
    /// panicking, and a free plain candidate is reused by later expansions.
    fn expansion_carrier(payload: AtomView<'_>, selected: AtomView<'_>) -> InferenceResult<Symbol> {
        let mut reserved = payload.get_all_symbols(true);
        reserved.extend(selected.get_all_symbols(true));
        let mut attempt = 0_u64;
        loop {
            if let Ok(symbol) = SymbolBuilder::new(NamespacedSymbol::parse(&format!(
                "idenso::expansion_selector_{attempt}"
            )))
            .with_attributes(Vec::<SymbolAttribute>::new())
            .build()
                && !reserved.contains(&symbol)
            {
                return Ok(symbol);
            }
            attempt = attempt.checked_add(1).ok_or_else(|| {
                TensorInferenceError::invalid("expansion selector symbol counter overflow")
            })?;
        }
    }

    /// Carry selection through arithmetic without opening original function
    /// leaves. The backend also treats those functions as opaque: a target in
    /// their arguments selects enclosing arithmetic, not the arguments themselves.
    /// Removing the carrier may open wrappers introduced by scope preparation;
    /// original functions contain no carrier and are reused without rebuilding.
    fn map_expansion_carrier(
        value: AtomView<'_>,
        target: AtomView<'_>,
        marker: Option<&Atom>,
        select_nested: bool,
    ) -> Option<Atom> {
        if let Some(marker) = marker {
            if matches!(value, AtomView::Fun(_)) {
                return (value == target || select_nested && value.contains(target))
                    .then(|| value.to_owned() * marker);
            }
        } else if value == target {
            return Some(Atom::one());
        }
        let children = |children: Vec<AtomView<'_>>| {
            let replaced = children
                .iter()
                .map(|&child| Self::map_expansion_carrier(child, target, marker, select_nested))
                .collect::<Vec<_>>();
            replaced.iter().any(Option::is_some).then(|| {
                children
                    .into_iter()
                    .zip(replaced)
                    .map(|(child, replaced)| replaced.unwrap_or_else(|| child.to_owned()))
                    .collect::<Vec<_>>()
            })
        };
        match value {
            AtomView::Add(sum) => {
                children(sum.iter().collect()).map(|terms| Atom::add_many(&terms))
            }
            AtomView::Mul(product) => {
                children(product.iter().collect()).map(|factors| Atom::mul_many(&factors))
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                children(vec![base, exponent]).map(|parts| parts[0].pow(&parts[1]))
            }
            AtomView::Fun(function) => children(function.iter().collect()).map(|arguments| {
                FunctionBuilder::new(function.get_symbol())
                    .add_args(&arguments)
                    .finish()
            }),
            _ => None,
        }
    }

    /// Whether `expression` has a negative rational power of a non-numeric
    /// base outside function arguments, which the polynomial backend leaves
    /// untouched.
    fn contains_symbolic_reciprocal(expression: AtomView<'_>) -> bool {
        let mut found = false;
        expression.visitor(&mut |value| {
            if let AtomView::Pow(power) = value {
                let (base, exponent) = power.get_base_exp();
                found |= !matches!(base, AtomView::Num(_))
                    && Rational::try_from(exponent)
                        .is_ok_and(|exponent| exponent.numerator().is_negative());
            }
            !found && !matches!(value, AtomView::Fun(_))
        });
        found
    }

    /// Expand scalar algebra while retaining the established tensor interface.
    /// A selected variable or polynomial mode uses the corresponding Symbolica
    /// operation. Pairs closed by a sum keep contracting their original factors:
    /// one which collides with a neighbouring factor, or with another term of
    /// an expanded power, first receives a fresh dummy name. An expanded power
    /// raises each closed term of its sum as one bracketed scalar, so the
    /// result keeps the multinomial size. Renamed terms can merge again after
    /// [`Self::canonize`]. A selected function also selects its renamed copies.
    /// Scopes which the selected expansion leaves alone stay unchanged. When a
    /// selected polynomial expansion changes the expression and its result
    /// divides by a non-numeric factor, the expansion is redone on the plain
    /// route, whose terms in the selector can be grouped differently. Debug
    /// builds check the user-normalizer contract.
    ///
    /// # Errors
    ///
    /// Rejects invalid source metadata, a non-indeterminate expansion target, a
    /// selector-symbol counter overflow, or an
    /// expansion that cannot retain unresolved tensor interfaces.
    ///
    /// # Example
    ///
    /// ```
    /// use idenso::tensor::SymbolicTensor;
    /// use idenso::tensor::inference::TensorInferenceError;
    /// use symbolica::atom::{Atom, AtomCore};
    /// use symbolica::symbol;
    ///
    /// let x = Atom::var(symbol!("docs::expansion_x"));
    /// let source = SymbolicTensor::infer((&x + 1).pow(2))?;
    /// let expanded = source.expanded(None, false)?;
    /// assert_eq!(expanded.expression(), &(&x * &x + Atom::num(2) * x + 1));
    /// # Ok::<(), TensorInferenceError>(())
    /// ```
    pub fn expanded(&self, var: Option<AtomView<'_>>, via_poly: bool) -> InferenceResult<Self> {
        self.explicit_index_occurrences()?;
        #[cfg(feature = "reference-cases")]
        let _phase = timing::scope(Phase::Materialization);
        if let Some(var) = var
            && !matches!(var, AtomView::Var(_) | AtomView::Fun(_))
        {
            return Err(TensorInferenceError::invalid(
                "Expansion must be done wrt an indeterminate",
            ));
        }
        if var.is_none()
            && !via_poly
            && !self.expression.as_view().needs_normalization()
            && self.expression.is_expanded::<Atom>(None)
        {
            return Ok(self.clone());
        }
        // Mark original function leaves before the single scoping pass. No
        // original function is opened or rebuilt by this marking operation.
        let (marked, marker) = if let Some(selected @ AtomView::Fun(_)) = var {
            let marker = Atom::var(Self::expansion_carrier(
                self.expression.as_view(),
                selected,
            )?);
            let marked = Self::map_expansion_carrier(
                self.expression.as_view(),
                selected,
                Some(&marker),
                !via_poly,
            );
            (marked, Some(marker))
        } else {
            (None, None)
        };
        let dummies = OnceCell::new();
        let separated = ExplicitIndexOccurrences::separate_expansion_scopes(
            marked.as_ref().unwrap_or(&self.expression).as_view(),
            marker.as_ref().map(Atom::as_view).or(var),
            via_poly,
            &mut || {
                dummies
                    .get_or_init(|| Self::reserved_dummies([self]))
                    .fresh_index()
            },
        );
        // If no scope changed, retain the backend's exact original input and
        // selector, including its usual coefficient factorization.
        let carrier = separated.is_some().then_some(marker).flatten();
        let source = separated.as_ref().unwrap_or(&self.expression);
        let expression = if let Some(target) = carrier.as_ref().map(Atom::as_view).or(var) {
            let expression = if via_poly {
                let expression = source.expand_via_poly::<i32, Atom>(Some(target.to_owned()));
                // Raising a polynomial with expression coefficients can divide
                // by its leading coefficient without cancelling it. Such a
                // reciprocal can share its base with a denominator of the
                // input, so any reciprocal in a changed result falls back to
                // the plain route, which never divides by a coefficient.
                if expression != *source && Self::contains_symbolic_reciprocal(expression.as_view())
                {
                    return self.expanded(var, false);
                }
                expression
            } else {
                source.expand_in(target)
            };
            match &carrier {
                Some(carrier) => {
                    let unmarked = Self::map_expansion_carrier(
                        expression.as_view(),
                        carrier.as_view(),
                        None,
                        false,
                    );
                    unmarked.unwrap_or(expression)
                }
                None => expression,
            }
        } else if via_poly {
            source.expand_via_poly::<i32, Atom>(None)
        } else {
            if self.structure.open_positions().is_empty()
                && let Some(expression) = SlotContraction::expand_scalar_sum(source.as_view())
            {
                if expression == self.expression {
                    return Ok(self.clone());
                }
                // Private intake proved every original leaf scalar and intrinsic,
                // before cancellation, and emitted only its normalized original.
                // This is the same successful scalar-algebra interface proof.
                return self.with_produced_expression(expression);
            }
            source.expand()
        };
        let mut result = self.expansion_result(self.borrowed_parts(), expression)?;
        // Cancellation can leave one term whose formerly closed pairs now
        // belong to the outermost product.
        if let Some(root) = ExplicitIndexOccurrences::root_product(result.expression.as_view())
            && let Some(occurrences) = result.proofs.occurrences.get_mut()
        {
            let occurrences = Arc::make_mut(occurrences);
            occurrences.branch_max(&root);
            occurrences.validate()?;
        }
        Ok(result)
    }

    /// Expansion distributes scalar algebra over the admitted leaves. Without
    /// open ports and with intrinsic normalization it keeps the established
    /// interface, even when the payload alone would not re-admit, for example
    /// a product of directed ports; other results take the algebra boundary.
    fn expansion_result(
        &self,
        parts: BorrowedTensorParts<'_>,
        expression: Atom,
    ) -> InferenceResult<Self> {
        if self.structure.open_positions().is_empty() && self.normalization_is_intrinsic() {
            let result = parts.with_established_expression(expression)?;
            TensorProofs::assert_normalization_preserves_interface(
                "expanded",
                &[&self.expression],
                &result.expression,
                Ok(&result.structure),
                None,
            );
            return Ok(result);
        }
        parts.with_algebra_result(expression)
    }

    /// Publish known scalar algebra which treats tensor leaves as opaque.
    /// Arbitrary user edits must use [`Self::with_checked_expression`]. Legal
    /// callbacks can rename dummies, so callback-bearing results account for
    /// emitted occurrences while transporting the authoritative free interface.
    pub fn with_algebra_result(&self, expression: Atom) -> InferenceResult<Self> {
        self.borrowed_parts().with_algebra_result(expression)
    }

    /// Publish a built-in result with this tensor's authoritative interface and
    /// construction history. The operation must preserve free-index identities,
    /// including encoded AUTO owners. Intrinsic operations must introduce no
    /// unreported dummy indices; callback-bearing outputs account newly written
    /// occurrences because legal callbacks may rename internal dummies.
    /// Changed payloads discard reduction caches; exact admitted region facts
    /// remain reusable for unchanged subexpressions.
    ///
    /// Arbitrary user transformations must use [`Self::with_checked_expression`]
    /// or [`Self::with_transformed_expression`] instead. Debug builds verify the
    /// producer contract without using observations as result metadata. See the
    /// [producer transport example](Self::with_produced_parts).
    ///
    /// # Errors
    /// Rejects an unpublished source without authoritative construction history.
    ///
    /// # Panics
    /// Panics in debug builds when the producer changes the admitted interface.
    pub fn with_produced_expression(&self, expression: Atom) -> InferenceResult<Self> {
        self.borrowed_parts().with_produced_expression(expression)
    }

    /// Publish a known interface transformation and apply its slot map to all
    /// retained construction history, including identities hidden by cancellation.
    /// The producer supplies the new logical port order and the same slot map
    /// used to rebuild the payload. Release builds do not infer its interface.
    ///
    /// # Errors
    /// Rejects unpublished metadata or a map producing invalid multiplicities.
    ///
    /// # Panics
    /// Panics in debug builds when the payload disagrees with the supplied interface.
    ///
    /// # Example
    ///
    /// An identity keeps admitted ports and history; a materializer additionally
    /// records a newly allocated internal contraction.
    ///
    /// ```
    /// use idenso::tensor::SymbolicTensor;
    /// use idenso::tensor::inference::TensorInferenceError;
    /// use spenso::structure::abstract_index::AbstractIndex;
    /// use spenso::structure::representation::{Euclidean, RepName};
    /// use spenso::structure::slot::IsAbstractSlot;
    /// use spenso::tensor_symbol;
    /// use symbolica::atom::Atom;
    /// use symbolica::function;
    /// let rep = Euclidean {}.new_rep(3).to_lib();
    /// let free = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(1));
    /// let dummy = rep.slot::<AbstractIndex, _>(AbstractIndex::Dummy(2));
    /// let left = tensor_symbol!("docs::producer_left");
    /// let right = tensor_symbol!("docs::producer_right");
    /// let source = SymbolicTensor::infer(function!(left, free.to_atom()))?;
    /// let weighted = source.with_produced_expression(Atom::num(2) * source.expression())?;
    /// let expanded = weighted.expression()
    ///     * function!(left, dummy.to_atom()) * function!(right, dummy.to_atom());
    /// let materialized = weighted.with_materialized_expression("known expansion", expanded)?;
    /// let transported = materialized.with_materialized_parts(
    ///     "known identity", materialized.expression().clone(),
    ///     materialized.structure().clone(), |slot| slot,
    /// )?;
    /// let zero = transported.with_produced_parts(
    ///     Atom::Zero, transported.structure().clone(), |slot| slot,
    /// )?;
    /// assert_eq!(zero.rank(), 1);
    /// assert!(zero.explicit_indices()?.any(|index| index == dummy.aind()));
    /// # Ok::<(), TensorInferenceError>(())
    /// ```
    pub fn with_produced_parts(
        &self,
        expression: Atom,
        structure: PartialStructure,
        mut map_slot: impl FnMut(Slot<LibraryRep, AbstractIndex>) -> Slot<LibraryRep, AbstractIndex>,
    ) -> InferenceResult<Self> {
        if expression != self.expression
            && !expression.is_zero()
            && !self.normalization_is_intrinsic()
        {
            return self.with_materialized_parts(
                "with_produced_parts",
                expression,
                structure,
                map_slot,
            );
        }
        let mut occurrences = self
            .explicit_index_occurrences()?
            .mapped(|slot| Ok::<_, TensorInferenceError>(map_slot(slot)))?;
        occurrences.branch_max(&ExplicitIndexOccurrences::from_interface(&structure));
        occurrences.validate()?;
        #[cfg(debug_assertions)]
        {
            #[cfg(test)]
            let _debug = DebugInference::enter();
            Self::assert_produced_interface("with_produced_parts", &expression, &structure);
            let observed = ExplicitIndexOccurrences::from_atom(
                expression.as_view(),
                &mut SlotMatcher::default(),
            );
            assert!(
                observed.is_bounded_by(&occurrences),
                "tensor producer introduced unreported index occurrences in with_produced_parts: {expression}"
            );
        }
        let mut result = Self::from_validated_parts(expression, structure, Arc::new(occurrences));
        result.proofs.inherit_interfaces([&self.proofs]);
        Ok(result)
    }

    #[cfg(debug_assertions)]
    fn assert_produced_interface(operation: &str, expression: &Atom, structure: &PartialStructure) {
        #[cfg(test)]
        let _debug = DebugInference::enter();
        let checked =
            InterfaceInference::merge_explicit_interface_sequence(slice::from_ref(structure))
                .and_then(|expected| {
                    Self::validate_observed_interface(
                        expression,
                        &expected,
                        LeafInference::for_declaration(&expected),
                    )
                });
        assert!(
            checked.is_ok(),
            "tensor producer contract violated in {operation}: result {expression}; admitted interface {structure:?}: {checked:?}"
        );
    }

    #[cfg(debug_assertions)]
    fn assert_constructor_interface(
        operation: &str,
        expression: &Atom,
        structure: &PartialStructure,
    ) {
        #[cfg(test)]
        let _debug = DebugInference::enter();
        // Constructors may declare a public axis order distinct from the
        // written order and may retain both ends of an internal contraction.
        // Compare the exposed incidences. Encoded AUTO indices still carry
        // literal owner and axis identities wherever the declaration does.
        let expected =
            InterfaceInference::merge_explicit_interface_sequence(slice::from_ref(structure));
        let checked = expected.and_then(|expected| {
            // Admission reads tensor-free input as a scalar without an
            // interface walk; its untagged calls are opaque data.
            if expression.is_zero()
                || (expected.logical_slots().is_empty()
                    && !InterfaceInference::contains_tensor_symbol(expression.as_view()))
            {
                return Ok(());
            }
            let observed = Self::observed_interface(
                expression,
                LeafInference::for_declaration(&expected),
                ObservationPurpose::Validate,
            )?;
            let mut ports = observed.interface.logical_slots();
            for port in expected.logical_slots() {
                let position = ports
                    .iter()
                    .position(|observed| {
                        port == *observed
                            || (port.rep() == observed.rep()
                                && matches!(
                                    (port.aind, observed.aind),
                                    (PartialIndex::Open(_), PartialIndex::Open(_))
                                ))
                    })
                    .ok_or_else(|| {
                        TensorInferenceError::invalid("result is missing a declared tensor port")
                    })?;
                ports.swap_remove(position);
            }
            if ports.is_empty() {
                Ok(())
            } else {
                Err(TensorInferenceError::invalid(
                    "result exposes undeclared tensor ports",
                ))
            }
        });
        assert!(
            checked.is_ok(),
            "tensor producer contract violated in {operation}: result {expression}; admitted interface {structure:?}: {checked:?}"
        );
    }

    /// Retain this interface after a built-in tensor identity. Normalization
    /// callbacks must preserve free indices; arbitrary user replacements use a
    /// checked boundary instead.
    pub fn with_rewritten_expression(&self, expression: Atom) -> InferenceResult<Self> {
        self.borrowed_parts().with_rewritten_expression(expression)
    }

    /// Infer an arbitrary transformation and report whether finishing contracted
    /// raw inferred ports. Presentation adapters use that fact to discard a
    /// stored-data descriptor, even when the input and final ranks are equal.
    pub fn with_transformed_expression(&self, expression: Atom) -> InferenceResult<(Self, bool)> {
        self.borrowed_parts()
            .with_transformed_expression(expression)
    }

    pub(crate) fn finish_mapped_result(
        &self,
        operation: &str,
        mut result: Self,
        mapped: TensorProofs,
    ) -> InferenceResult<Self> {
        result.proofs.inherit_interfaces([&mapped]);
        if result.expression != self.expression && !self.normalization_is_intrinsic() {
            let mut occurrences = InterfaceInference::produced_occurrences(
                result.expression.as_view(),
                &mut SlotMatcher::default(),
                &result.proofs,
            )?;
            occurrences.branch_max(result.explicit_index_occurrences()?);
            occurrences.validate()?;
            result.proofs.occurrences = OnceLock::from(Arc::new(occurrences));
        }
        TensorProofs::assert_normalization_preserves_interface(
            operation,
            &[&self.expression],
            &result.expression,
            Ok(&result.structure),
            None,
        );
        Ok(result)
    }

    /// Change four-dimensional Lorentz ports in the expression and its logical interface.
    /// The same map transports retained index history, including identities hidden
    /// by cancellation. Rebuilding callbacks must preserve the admitted tensor meaning.
    /// See the [metadata transport example](Self::with_renamed_indices).
    ///
    /// # Errors
    ///
    /// Returns the dimension-encoding, admission, and multiplicity errors from
    /// [`Self::map_representation_dimensions`].
    pub fn with_lorentz_dimension(&self, dimension: Dimension) -> InferenceResult<Self> {
        self.map_representation_dimensions(|representation, current| {
            if representation == (Minkowski {}).into() && current == Atom::num(4).as_view() {
                dimension.to_symbolic()
            } else {
                current.to_owned()
            }
        })
    }

    /// Map representation dimensions in both the payload and its logical ports.
    /// Portable encoded dimensions are presented to the mapper in source notation;
    /// only the resulting dimension is encoded again. Index identities are untouched.
    /// Each distinct representation/dimension pair is mapped once, including ports
    /// carried only by a ranked zero. Scalar coefficients outside representations
    /// are not substituted by this operation.
    /// See the [example](Self::with_renamed_indices).
    ///
    /// # Errors
    ///
    /// Rejects invalid source metadata, dimensions that cannot be encoded reversibly,
    /// invalid mapped slots, or newly incompatible index multiplicities.
    pub fn map_representation_dimensions(
        &self,
        mut map: impl FnMut(LibraryRep, AtomView<'_>) -> Atom,
    ) -> InferenceResult<Self> {
        self.explicit_index_occurrences()?;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(false, true);
        let mut dimensions = HashMap::<Representation<LibraryRep>, Dimension>::new();
        let mapped_representation = RefCell::new(|representation: Representation<LibraryRep>| {
            let dimension = if let Some(dimension) = dimensions.get(&representation) {
                *dimension
            } else {
                let encoded = representation.dim.to_symbolic();
                let source = encoded
                    .get_symbol()
                    .and_then(CookSettings::dimension_payload)
                    .unwrap_or(encoded);
                let mapped = map(representation.rep, source.as_view());
                let compact = representation.rep.to_symbolic([mapped]);
                let compact = cooking
                    .try_cook(compact.as_view())
                    .map_err(TensorInferenceError::MappedDimensionEncoding)?;
                let dimension = Representation::<LibraryRep>::try_from(compact.as_view())
                    .map_err(|error| TensorInferenceError::Slot {
                        operation: "invalid mapped dimension",
                        source: error,
                    })?
                    .dim;
                dimensions.insert(representation, dimension);
                dimension
            };
            Ok::<_, TensorInferenceError>(Representation {
                rep: representation.rep,
                dim: dimension,
            })
        });
        let mapped_structure = |structure: &PartialStructure| {
            structure
                .logical_slots()
                .into_iter()
                .map(|slot| Ok((mapped_representation.borrow_mut())(slot.rep())?.slot(slot.aind)))
                .collect::<InferenceResult<Vec<_>>>()
                .map(PartialStructure::from_logical_slots)
        };
        let structure = mapped_structure(&self.structure)?;
        let mut slots = SlotMatcher::default();
        let (expression, mapped) = self.map_admitted_expression(
            |view| {
                let (representation, index) = if let SlotMatch::Explicit(slot) =
                    slots.classify(view)
                {
                    let representation =
                        Representation::<LibraryRep>::try_from(view).map_err(|source| {
                            TensorInferenceError::Slot {
                                operation: "invalid source representation",
                                source,
                            }
                        })?;
                    (representation, Some(slot.index().to_owned()))
                } else if let Ok(representation) = slots.parse_representation::<LibraryRep>(view) {
                    (representation, None)
                } else {
                    return Ok(None);
                };
                let mapped = (mapped_representation.borrow_mut())(representation)?;
                Ok(Some(if mapped != representation {
                    mapped.to_symbolic(index)
                } else {
                    view.to_owned()
                }))
            },
            mapped_structure,
        )?;
        let source_occurrences = self.explicit_index_occurrences()?;
        let occurrences = source_occurrences.mapped(|slot| {
            Ok::<_, TensorInferenceError>(
                (mapped_representation.borrow_mut())(slot.rep())?.slot(slot.aind()),
            )
        })?;
        occurrences.validate()?;
        let metadata_unchanged = occurrences
            .oriented_entries()
            .eq(source_occurrences.oriented_entries());
        if expression == self.expression && structure == self.structure && metadata_unchanged {
            // An unchanged carrier retains its established proofs while
            // still validating unchecked inputs and normalizing unfinished parts.
            return self.clone().finish_parts();
        }
        let structure = InterfaceInference::merge_explicit_interface_sequence(&[structure])?
            .canonicalize_open_ports();
        let mut value = Self::from_validated_parts(expression, structure, Arc::new(occurrences));
        if value.expression == self.expression {
            value.is_metric = self.is_metric;
            value.is_composite = self.is_composite;
        }
        self.finish_mapped_result("map_representation_dimensions", value, mapped)
    }

    /// Scope explicit indices and their logical ports; AUTO identities are unchanged.
    /// With `dummies_only`, protect the established external slots from scoping;
    /// a sum which closes such a name internally scopes that pair together.
    /// See the [example](Self::with_renamed_indices).
    ///
    /// # Errors
    ///
    /// Rejects invalid source metadata or invalid resulting explicit-index multiplicities.
    pub fn wrap_indices(&self, scope: Symbol, dummies_only: bool) -> InferenceResult<Self> {
        self.explicit_index_occurrences()?;
        let mut structure = self.structure.logical_slots();
        let external = structure
            .iter()
            .filter_map(|slot| match slot.aind {
                PartialIndex::Explicit(index) => Some(slot.rep().slot(index)),
                PartialIndex::Open(_) => None,
            })
            .collect::<HashSet<_>>();
        let mut matcher = SlotMatcher::default();
        let map_slot = |slot: Slot<LibraryRep, AbstractIndex>| {
            if !matches!(slot.aind(), AbstractIndex::Open { .. })
                && (!dummies_only || !external.contains(&slot))
            {
                slot.rep().slot(slot.aind().scoped(scope))
            } else {
                slot
            }
        };
        let (expression, mapped) = self.map_admitted_expression(
            |view| {
                if view.get_symbol().is_some_and(|symbol| symbol.is_scalar()) {
                    return Ok(Some(view.to_owned()));
                }
                Ok(matcher
                    .parse::<LibraryRep, AbstractIndex>(view)
                    .ok()
                    .map(|slot| {
                        let mapped = map_slot(slot);
                        if mapped == slot {
                            view.to_owned()
                        } else {
                            mapped.to_atom()
                        }
                    }))
            },
            |interface| {
                Ok(PartialStructure::from_logical_slots(
                    interface
                        .logical_slots()
                        .into_iter()
                        .map(|slot| match slot.aind {
                            PartialIndex::Explicit(index) => {
                                let mapped = map_slot(slot.rep().slot(index));
                                mapped.rep().slot(PartialIndex::Explicit(mapped.aind()))
                            }
                            PartialIndex::Open(_) => slot,
                        }),
                ))
            },
        )?;
        if !dummies_only {
            for slot in &mut structure {
                if let PartialIndex::Explicit(index) = slot.aind {
                    slot.aind = PartialIndex::Explicit(index.scoped(scope));
                }
            }
        }
        let structure = PartialStructure::from_logical_slots(structure);
        let occurrences = self.explicit_index_occurrences()?.mapped(|slot| {
            let protected = external.contains(&slot);
            Ok::<_, TensorInferenceError>(if !dummies_only || !protected {
                slot.rep().slot(slot.aind().scoped(scope))
            } else {
                slot
            })
        })?;
        occurrences.validate()?;
        let mut result = Self::from_validated_parts(expression, structure, Arc::new(occurrences));
        result.proofs.normalized_product = self.proofs.normalized_product.clone();
        self.finish_mapped_result("wrap_indices", result, mapped)
    }

    /// Carry the protected external boundary through sums and closed powers.
    /// Only changed ancestors are reconstructed, so discovery never replays a
    /// user normalizer before the actual slot transformation. Only arguments
    /// which inference reads as indices are scoped: an untagged function is an
    /// opaque scalar whose arguments are data, even when they look like slots.
    fn wrap_unprotected_indices(
        expression: AtomView<'_>,
        scope: Symbol,
        protected: &HashSet<Slot<LibraryRep, AbstractIndex>>,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        if expression.get_symbol().is_some_and(|head| head.is_scalar()) {
            return None;
        }
        if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(expression) {
            if protected.contains(&slot) {
                return None;
            }
            let wrapped = AbstractIndex::wrap_expression(expression, scope, |_| true);
            return (wrapped.as_view() != expression).then_some(wrapped);
        }
        let boundary = (!protected.is_empty()
            && (matches!(expression, AtomView::Add(_))
                || Self::is_closed_scalar_power(expression, slots)))
        .then(|| {
            let exported = ExplicitIndexOccurrences::from_atom(expression, slots);
            protected
                .iter()
                .filter(|slot| exported.count(**slot) != 0)
                .copied()
                .collect::<HashSet<_>>()
        });
        let protected = boundary.as_ref().unwrap_or(protected);
        let mut children =
            |children: Vec<AtomView<'_>>,
             indexed: &dyn Fn(AtomView<'_>, &mut SlotMatcher) -> bool| {
                let replaced = children
                    .iter()
                    .map(|&child| {
                        indexed(child, slots)
                            .then(|| Self::wrap_unprotected_indices(child, scope, protected, slots))
                            .flatten()
                    })
                    .collect::<Vec<_>>();
                replaced.iter().any(Option::is_some).then(|| {
                    children
                        .into_iter()
                        .zip(replaced)
                        .map(|(child, replaced)| replaced.unwrap_or_else(|| child.to_owned()))
                        .collect::<Vec<_>>()
                })
            };
        let all = |_: AtomView<'_>, _: &mut SlotMatcher| true;
        match expression {
            AtomView::Add(sum) => children(sum.iter().collect(), &all).map(Atom::add_many),
            AtomView::Mul(product) => children(product.iter().collect(), &all).map(Atom::mul_many),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                children(vec![base, exponent], &all).map(|parts| parts[0].pow(&parts[1]))
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                children(function.iter().collect(), &|argument, slots| {
                    InterfaceInference::explicit_index_argument(symbol, argument, slots)
                })
                .map(|arguments| FunctionBuilder::new(symbol).add_args(&arguments).finish())
            }
            _ => None,
        }
    }

    /// Whether index observation closes this power's base, as it does for a
    /// scalar `bracket(…)^k`: the base's observed interface is scalar, so the
    /// base owns its dummy pairs for every exponent. A tensor-valued base
    /// instead repeats its ports into the surrounding product. Observation
    /// neither rebuilds the base nor invokes its normalization callbacks.
    fn is_closed_scalar_power(node: AtomView<'_>, slots: &mut SlotMatcher) -> bool {
        let AtomView::Pow(power) = node else {
            return false;
        };
        let (base, _) = power.get_base_exp();
        let mut inference = InterfaceInference {
            slots: mem::take(slots),
            leaf_inference: LeafInference::Observe,
            ..InterfaceInference::default()
        };
        let interface = inference.infer_view(base).and_then(|interface| {
            InterfaceInference::merge_explicit_interface_sequence(&[interface])
        });
        *slots = inference.slots;
        interface.is_ok_and(|interface| interface.canonical().is_scalar())
    }

    /// Cook both encoded indices and their logical interface, contracting collisions.
    /// See the [example](Self::with_renamed_indices).
    ///
    /// # Errors
    ///
    /// Rejects invalid source metadata, unsupported cooking payloads, malformed cooked
    /// slots, or resulting explicit-index multiplicities greater than two.
    pub fn with_cooked_indices(&self, settings: &CookSettings) -> InferenceResult<Self> {
        self.explicit_index_occurrences()?;
        let cook = |atom: &Atom| {
            settings
                .try_cook_indices(atom.as_view())
                .map_err(TensorInferenceError::Cooking)
        };
        let map_interface = |interface: &PartialStructure| {
            let slots = interface
                .logical_slots()
                .into_iter()
                .map(|slot| {
                    if matches!(slot.aind, PartialIndex::Open(_)) {
                        return Ok(slot);
                    }
                    let cooked = cook(&composition::port_atom(slot))?;
                    let cooked = Slot::<LibraryRep, AbstractIndex>::try_from(cooked.as_view())
                        .map_err(|error| TensorInferenceError::Slot {
                            operation: "cannot cook interface slot",
                            source: error,
                        })?;
                    Ok(cooked.rep().slot(PartialIndex::Explicit(cooked.aind())))
                })
                .collect::<InferenceResult<Vec<_>>>()?;
            InterfaceInference::merge_explicit_interface_sequence(&[
                PartialStructure::from_logical_slots(slots),
            ])
            .map(|structure| structure.canonicalize_open_ports())
        };
        let structure = map_interface(&self.structure)?;
        let (expression, mapped) = self.map_admitted_expression(
            |view| {
                let AtomView::Fun(function) = view else {
                    return Ok(None);
                };
                if function.get_symbol().is_scalar() {
                    return Ok(Some(view.to_owned()));
                }
                if !function.get_symbol().has_tag(&SPENSO_TAG.representation)
                    || !(1..=2).contains(&function.get_nargs())
                {
                    return Ok(None);
                }
                let cooked = settings
                    .try_cook_indices(view)
                    .map_err(TensorInferenceError::Cooking)?;
                Ok((cooked.as_view() != view).then_some(cooked))
            },
            map_interface,
        )?;
        let occurrences = self.explicit_index_occurrences()?.mapped(|slot| {
            let cooked = cook(&slot.to_atom())?;
            Slot::<LibraryRep, AbstractIndex>::try_from(cooked.as_view()).map_err(|error| {
                TensorInferenceError::Slot {
                    operation: "cannot cook retained index",
                    source: error,
                }
            })
        })?;
        occurrences.validate()?;
        let value = Self::from_validated_parts(expression, structure, Arc::new(occurrences));
        self.finish_mapped_result("with_cooked_indices", value, mapped)
    }

    /// Rename symbols in the payload and its interface. The renaming must be
    /// injective and avoid symbols already present, so the established
    /// interface and validation carry over unchanged.
    pub(crate) fn with_renamed_symbols(
        &self,
        renaming: &HashMap<Symbol, Symbol>,
    ) -> InferenceResult<Self> {
        self.explicit_index_occurrences()?;
        let rename = |atom: AtomView<'_>| {
            atom.replace_map(|node, _context, out| {
                if let AtomView::Var(variable) = node
                    && let Some(target) = renaming.get(&variable.get_symbol())
                {
                    **out = Atom::var(*target);
                }
            })
        };
        let map_interface = |interface: &PartialStructure| {
            let slots = interface
                .logical_slots()
                .into_iter()
                .map(|slot| {
                    let renamed = rename(composition::port_atom(slot).as_view());
                    if matches!(slot.aind, PartialIndex::Open(_)) {
                        return Representation::<LibraryRep>::try_from(renamed.as_view())
                            .map(|representation| representation.slot(slot.aind))
                            .map_err(|error| TensorInferenceError::Slot {
                                operation: "cannot rename interface representation",
                                source: error,
                            });
                    }
                    let renamed = Slot::<LibraryRep, AbstractIndex>::try_from(renamed.as_view())
                        .map_err(|error| TensorInferenceError::Slot {
                            operation: "cannot rename interface slot",
                            source: error,
                        })?;
                    Ok(renamed.rep().slot(PartialIndex::Explicit(renamed.aind())))
                })
                .collect::<InferenceResult<Vec<_>>>()?;
            Ok(PartialStructure::from_logical_slots(slots))
        };
        let structure = map_interface(&self.structure)?;
        let (expression, mapped) = self.map_admitted_expression(
            |view| {
                Ok(if let AtomView::Var(variable) = view {
                    renaming
                        .get(&variable.get_symbol())
                        .map(|target| Atom::var(*target))
                } else {
                    None
                })
            },
            map_interface,
        )?;
        let occurrences = self.explicit_index_occurrences()?.mapped(|slot| {
            let renamed = rename(slot.to_atom().as_view());
            Slot::<LibraryRep, AbstractIndex>::try_from(renamed.as_view()).map_err(|error| {
                TensorInferenceError::Slot {
                    operation: "cannot rename retained index",
                    source: error,
                }
            })
        })?;
        occurrences.validate()?;
        let mut result = Self::from_validated_parts(expression, structure, Arc::new(occurrences));
        result.proofs.normalized_product = self.proofs.normalized_product.clone();
        self.finish_mapped_result("with_renamed_symbols", result, mapped)
    }

    /// Carry admitted metadata through a producer-certified index relabeling.
    ///
    /// `expression` must implement this simultaneous relabeling and may bind
    /// scalar owner/momentum arguments without changing tensor signatures. The
    /// producer must not introduce other indices, change dimensions, or run an
    /// interface-changing callback. Unmapped indices retain their identity.
    /// Injectivity includes retained identities absent from the simplified
    /// payload. Pure relabeling uses [`Self::renamed_indices`]; additional
    /// user-supplied changes are checked at this boundary.
    ///
    /// # Errors
    ///
    /// Rejects invalid source metadata, noninjective maps (including collisions with
    /// retained hidden identities), unresolved target indices, or excessive multiplicity.
    ///
    /// # Example
    ///
    /// The producer supplies the payload with the same index substitution as the map.
    ///
    /// ```
    /// use std::collections::HashMap;
    /// use idenso::CookSettings;
    /// use idenso::tensor::SymbolicTensor;
    /// use idenso::tensor::inference::TensorInferenceError;
    /// use spenso::structure::abstract_index::AbstractIndex;
    /// use spenso::structure::dimension::Dimension;
    /// use spenso::structure::representation::{Minkowski, RepName};
    /// use spenso::structure::slot::IsAbstractSlot;
    /// use spenso::tensor_symbol;
    /// use symbolica::atom::Atom;
    /// use symbolica::{function, symbol};
    ///
    /// let rep = Minkowski {}.new_rep(4).to_lib();
    /// let old = AbstractIndex::Symbol(symbol!("docs::old_index").into());
    /// let new = AbstractIndex::Symbol(symbol!("docs::new_index").into());
    /// let head = tensor_symbol!("docs::renamed_vector");
    /// let old_port = rep.slot::<AbstractIndex, _>(old).to_atom();
    /// let new_port = rep.slot::<AbstractIndex, _>(new).to_atom();
    /// let source = SymbolicTensor::infer(function!(head, old_port))?;
    /// let payload = function!(head, new_port);
    /// let renamed = source.with_renamed_indices(payload, &HashMap::from([(old, new)]))?;
    /// assert_eq!(renamed.rank(), 1);
    /// let extended = renamed.with_lorentz_dimension(Dimension::Concrete(6))?;
    /// let mapped = renamed.map_representation_dimensions(|_, _| Atom::num(6))?;
    /// assert!(mapped.same_value(&extended));
    /// assert_eq!(renamed.wrap_indices(symbol!("docs::scope"), false)?.rank(), 1);
    /// assert!(renamed.with_cooked_indices(&CookSettings::indices())?.same_value(&renamed));
    /// let relabelled = source.rename_symbols(&HashMap::from([
    ///     (symbol!("docs::old_index"), symbol!("docs::new_index")),
    /// ]))?;
    /// assert!(relabelled.same_value(&renamed));
    /// # Ok::<(), TensorInferenceError>(())
    /// ```
    pub fn with_renamed_indices(
        &self,
        expression: Atom,
        renaming: &HashMap<AbstractIndex, AbstractIndex>,
    ) -> InferenceResult<Self> {
        let renamed = self.renamed_indices(renaming)?;
        if expression == renamed.expression {
            return Ok(renamed);
        }
        // The supplied payload may contain arbitrary scalar substitutions. Its
        // changed leaves have no producer provenance, so this is a user boundary.
        renamed.with_checked_expression(expression)
    }

    /// Simultaneously rename admitted indices, carrying every regional and raw
    /// leaf interface through the same construction. Hidden index history is
    /// mapped too. See [`Self::with_renamed_indices`] for a complete example.
    ///
    /// # Errors
    /// Rejects noninjective maps, collisions, and unresolved target indices.
    pub fn renamed_indices(
        &self,
        renaming: &HashMap<AbstractIndex, AbstractIndex>,
    ) -> InferenceResult<Self> {
        let occurrences = self.renamed_occurrences(renaming)?;
        let structure = Self::renamed_interface(&self.structure, renaming);
        let (expression, mapped) = self.map_admitted_indices(
            |view| {
                if view.get_symbol().is_some_and(|head| head.is_scalar()) {
                    return Ok(Some(view.to_owned()));
                }
                Ok(Slot::<LibraryRep, AbstractIndex>::try_from(view)
                    .ok()
                    .and_then(|slot| {
                        renaming
                            .get(&slot.aind())
                            .map(|target| slot.rep().slot::<AbstractIndex, _>(*target).to_atom())
                    }))
            },
            |interface| Ok(Self::renamed_interface(interface, renaming)),
        )?;
        let mut result = Self::from_validated_parts(expression, structure, Arc::new(occurrences));
        result.proofs.normalized_product = self.proofs.normalized_product.clone();
        self.finish_mapped_result("renamed_indices", result, mapped)
    }

    /// Map an established producer's expression and its admitted interfaces together.
    ///
    /// The mapper runs once per source node, after its children have been mapped.
    /// It must implement `renaming` and may bind scalar metadata or replace a
    /// tensor head by a producer-known head with the same port signature.
    /// A rank-one leaf may become a scalar linear combination of rank-one
    /// constructors carrying its same mapped raw port. Other tensor slot roles
    /// must remain unchanged. It must not introduce new index identities or
    /// change dimensions. Normalization callbacks must preserve the mapped free
    /// interface. Unknown user transformations use [`Self::with_checked_expression`].
    /// Hidden index history and logical port order follow `renaming`, including zero.
    ///
    /// # Errors
    /// Returns mapping errors and rejects noninjective or unresolved index maps.
    ///
    /// # Example
    /// ```
    /// use std::collections::HashMap;
    /// use idenso::tensor::SymbolicTensor;
    /// use idenso::tensor::inference::TensorInferenceError;
    /// use symbolica::atom::Atom;
    /// let source = SymbolicTensor::infer(Atom::num(3))?;
    /// let mapped = source.map_produced_expression(&HashMap::new(), |_| {
    ///     Ok::<_, TensorInferenceError>(None)
    /// })?;
    /// assert!(mapped.same_value(&source));
    /// # Ok::<(), TensorInferenceError>(())
    /// ```
    pub fn map_produced_expression<E>(
        &self,
        renaming: &HashMap<AbstractIndex, AbstractIndex>,
        node_map: impl for<'a> FnMut(AtomView<'a>) -> Result<Option<Atom>, E>,
    ) -> Result<Self, E>
    where
        E: From<TensorInferenceError>,
    {
        let occurrences = self.renamed_occurrences(renaming)?;
        let structure = Self::renamed_interface(&self.structure, renaming);
        let (expression, mapped) = self
            .map_admitted_expression_bottom_up(node_map, |interface| {
                Ok(Self::renamed_interface(interface, renaming))
            })?;
        let mut result = Self::from_validated_parts(expression, structure, Arc::new(occurrences));
        result.proofs.normalized_product = self.proofs.normalized_product.clone();
        self.finish_mapped_result("map_produced_expression", result, mapped)
            .map_err(E::from)
    }

    fn renamed_occurrences(
        &self,
        renaming: &HashMap<AbstractIndex, AbstractIndex>,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        let source = self.explicit_index_occurrences()?;
        let mut indices = source
            .entries()
            .map(|(slot, _)| slot.aind())
            .collect::<HashSet<_>>();
        indices.extend(self.structure.logical_slots().into_iter().filter_map(
            |slot| match slot.aind {
                PartialIndex::Explicit(index) => Some(index),
                PartialIndex::Open(_) => None,
            },
        ));
        indices.extend(renaming.keys().copied());
        let mut targets = HashSet::with_capacity(indices.len());
        for index in indices {
            let target = renaming.get(&index).copied().unwrap_or(index);
            if !targets.insert(target) {
                return Err(TensorInferenceError::invalid(
                    "index renaming must be injective and avoid existing identities",
                ));
            }
            if CookSettings::index_has_open_port(target)? {
                return Err(TensorInferenceError::invalid(
                    "index renaming cannot introduce unresolved ports",
                ));
            }
        }
        let occurrences = source.mapped(|slot| {
            Ok::<_, TensorInferenceError>(
                slot.rep()
                    .slot(renaming.get(&slot.aind()).copied().unwrap_or(slot.aind())),
            )
        })?;
        occurrences.validate()?;
        Ok(occurrences)
    }

    fn renamed_interface(
        interface: &PartialStructure,
        renaming: &HashMap<AbstractIndex, AbstractIndex>,
    ) -> PartialStructure {
        PartialStructure::from_logical_slots(interface.logical_slots().into_iter().map(
            |mut slot| {
                if let PartialIndex::Explicit(index) = slot.aind {
                    slot.aind =
                        PartialIndex::Explicit(renaming.get(&index).copied().unwrap_or(index));
                }
                slot
            },
        ))
    }

    /// Simultaneously rename symbols without merging distinct source identities.
    /// Swaps are allowed; an unmapped existing target or a noninjective map is not.
    /// Logical ports, including those of zero, follow the same substitution.
    /// Normalizers triggered by rebuilding the payload must preserve that mapped
    /// interface; observing their result does not replay them.
    /// See the [example](Self::with_renamed_indices).
    ///
    /// # Errors
    ///
    /// Rejects noninjective renaming, collisions with existing symbols (including hidden
    /// index history), invalid source metadata, or malformed mapped slots.
    pub fn rename_symbols(&self, renaming: &HashMap<Symbol, Symbol>) -> InferenceResult<Self> {
        self.explicit_index_occurrences()?;
        let mut present = self.expression.get_all_symbols(true);
        for slot in self.structure.logical_slots() {
            present.extend(composition::port_atom(slot).get_all_symbols(true));
        }
        for (slot, _) in self.explicit_index_occurrences()?.entries() {
            present.extend(slot.to_atom().get_all_symbols(true));
        }
        let mut targets = HashSet::new();
        for (source, target) in renaming {
            if !targets.insert(*target) {
                return Err(TensorInferenceError::invalid(
                    "symbol renaming must be injective",
                ));
            }
            if source != target && present.contains(target) && !renaming.contains_key(target) {
                return Err(TensorInferenceError::invalid(
                    "symbol renaming targets an existing unmapped symbol",
                ));
            }
        }
        let mut value = self.with_renamed_symbols(renaming)?;
        if value.expression == self.expression {
            value.is_metric = self.is_metric;
            value.is_composite = self.is_composite;
        }
        Ok(value)
    }

    /// Multiply by an independently scoped scalar expression.
    /// Dummy indices of the scalar are renamed when they collide with this tensor.
    /// Unlike [`ScalarMul::scalar_mul`](spenso::algebra::ScalarMul::scalar_mul),
    /// this fallible entry point reports invalid scalar input with a typed error.
    ///
    /// # Errors
    ///
    /// Rejects invalid source metadata, a non-scalar factor, or invalid scalar syntax.
    /// See the [arithmetic example](SymbolicTensor).
    pub fn try_scalar_mul(&self, rhs: &Atom) -> InferenceResult<Self> {
        let scalar = Self::infer(rhs.clone())?;
        if !scalar.is_scalar() {
            return Err(TensorInferenceError::invalid(
                "scalar multiplication requires a scalar factor",
            ));
        }
        let source = ExplicitIndexOccurrences::product(&[self])?;
        let dummies = Self::reserved_dummies([self, &scalar]);
        let renaming = scalar
            .explicit_index_occurrences()?
            .entries()
            .filter(|(slot, count)| *count == 2 && source.count(*slot) != 0)
            .map(|(slot, _)| (slot.aind(), dummies.fresh_index()))
            .collect::<HashMap<_, _>>();
        if renaming.is_empty() {
            return Ok(self.multiply(&scalar)?);
        }
        let scalar = scalar.renamed_indices(&renaming)?;
        Ok(self.multiply(&scalar)?)
    }

    /// Add tensors with compatible interfaces. Named results use canonical axes;
    /// unresolved occurrences retain their established logical correspondence.
    /// See the [example](SymbolicTensor).
    ///
    /// # Errors
    ///
    /// Rejects invalid operand metadata or incompatible logical interfaces.
    pub fn try_add(&self, right: &Self) -> InferenceResult<Self> {
        self.explicit_index_occurrences()?;
        right.explicit_index_occurrences()?;
        if !InterfaceInference::additive_interfaces_match(&self.structure, &right.structure) {
            return Err(TensorInferenceError::invalid(
                "addition requires compatible tensor interfaces",
            ));
        }
        self.additive_result(right, &self.expression + &right.expression)
    }

    /// Publish a normalized sum of operands with matching interfaces.
    fn additive_result(&self, right: &Self, expression: Atom) -> InferenceResult<Self> {
        // Each summand binds its internal pairs. Only the shared public ports
        // remain visible to subsequent products, including a ranked zero,
        // unless normalization collapsed the sum into one product.
        let mut occurrences = self.explicit_index_occurrences()?.scoped(&self.structure);
        occurrences.branch_max(&right.explicit_index_occurrences()?.scoped(&right.structure));
        if let Some(root) = ExplicitIndexOccurrences::root_product(expression.as_view()) {
            occurrences.branch_max(&root);
            occurrences.validate()?;
        }
        let mut result = Self::from_validated_parts(
            expression,
            self.structure.clone().into_algebra_result(),
            Arc::new(occurrences),
        );
        result
            .proofs
            .inherit_interfaces([&self.proofs, &right.proofs]);
        let _ = result.proofs.normalized_product.set(
            self.borrowed_parts().product_is_normalized()
                && right.borrowed_parts().product_is_normalized(),
        );
        Ok(result)
    }

    /// Subtract tensors with compatible logical interfaces, retaining typed zero.
    /// See the [example](SymbolicTensor).
    ///
    /// # Errors
    ///
    /// Rejects invalid operand metadata or incompatible logical interfaces.
    pub fn try_sub(&self, right: &Self) -> InferenceResult<Self> {
        self.explicit_index_occurrences()?;
        right.explicit_index_occurrences()?;
        if !InterfaceInference::additive_interfaces_match(&self.structure, &right.structure) {
            return Err(TensorInferenceError::invalid(
                "subtraction requires compatible tensor interfaces",
            ));
        }
        self.additive_result(right, &self.expression - &right.expression)
    }

    /// Validate index multiplicity and the scopes of chain and trace placeholders.
    pub fn validate_atom(atom: &Atom) -> InferenceResult<()> {
        InterfaceInference::scoped_occurrences(atom.as_view(), &mut SlotMatcher::default())?
            .validate()?;
        Ok(())
    }

    /// Check a rewrite against an established interface, retaining typed zeros.
    ///
    /// Index substitution can invoke a user normalizer which changes tensor rank;
    /// inspect the resulting syntax without invoking that normalizer again with
    /// synthetic indices. Constructors still materialize compact callback leaves.
    pub fn validate_interface(atom: &Atom, structure: &PartialStructure) -> InferenceResult<()> {
        Self::validate_observed_interface(atom, structure, LeafInference::Observe).map(|_| ())
    }

    /// Check network storage indices without discarding encoded open-port owners.
    #[cfg(test)]
    pub(crate) fn validate_encoded_interface(
        atom: &Atom,
        structure: &PartialStructure,
    ) -> InferenceResult<()> {
        Self::validate_observed_interface(atom, structure, LeafInference::ObserveEncoded)
            .map(|_| ())
    }

    pub(super) fn validate_observed_interface(
        atom: &Atom,
        structure: &PartialStructure,
        policy: LeafInference,
    ) -> InferenceResult<(PartialStructure, Option<ExplicitIndexOccurrences>)> {
        Self::validate_observed_parts(atom, structure, policy, ObservationPurpose::Validate)
            .map(|observed| (observed.interface, observed.occurrences))
    }

    fn validate_observed_parts(
        atom: &Atom,
        structure: &PartialStructure,
        policy: LeafInference,
        purpose: ObservationPurpose,
    ) -> InferenceResult<ObservedTensor> {
        if atom.as_view().is_zero() {
            return Ok(ObservedTensor {
                interface: structure.clone(),
                occurrences: None,
                interfaces: AdmittedInterfaces::default(),
                leaf_interfaces: Default::default(),
            });
        }
        let observed = match Self::observed_interface(atom, policy, purpose) {
            Err(TensorInferenceError::Invalid(_))
                if structure.canonical().order() != 0
                    && !InterfaceInference::has_structured_syntax(atom.as_view()) =>
            {
                // A nonzero scalar cannot replace an established tensor. Keep
                // that interface error ahead of malformed scalar scope errors;
                // scalar destinations still require full syntax admission.
                ObservedTensor {
                    interface: PartialStructure::from_logical_slots([]),
                    occurrences: None,
                    interfaces: AdmittedInterfaces::default(),
                    leaf_interfaces: Default::default(),
                }
            }
            result => result?,
        };
        let inferred = &observed.interface;
        if !InterfaceInference::additive_interfaces_match(structure, inferred)
            && !(matches!(
                policy,
                LeafInference::ObserveEncoded | LeafInference::ObserveStorage
            ) && InterfaceInference::unambiguous_anonymous_ports(structure)
                && InterfaceInference::unambiguous_anonymous_ports(inferred)
                && structure.canonical() == inferred.canonical())
        {
            return Err(TensorInferenceError::invalid(
                "transformed expression does not preserve a compatible tensor interface",
            ));
        }
        Ok(observed)
    }

    /// Observe encoded port order without materializing compact callback leaves.
    fn observed_interface(
        atom: &Atom,
        policy: LeafInference,
        purpose: ObservationPurpose,
    ) -> InferenceResult<ObservedTensor> {
        // Scalar results share syntax admission with constructor results. In
        // particular an untagged structural carrier cannot bypass validation
        // merely because the graph's Tagged policy treats it as opaque.
        let mut inference = InterfaceInference {
            leaf_inference: policy,
            purpose,
            ..InterfaceInference::default()
        };
        let (interface, occurrences) = inference.infer_observed(atom.as_view())?;
        let interface = if policy == LeafInference::ObserveStorage {
            interface
        } else {
            InterfaceInference::merge_explicit_interface_sequence(&[interface])?
        };
        Ok(ObservedTensor {
            interface: interface.canonicalize_open_ports(),
            occurrences: Some(occurrences),
            interfaces: inference.admitted_interfaces,
            leaf_interfaces: Arc::new(Mutex::new(inference.admitted_leaf_interfaces)),
        })
    }
}

impl InterfaceInference {
    // These structural builtins cannot invoke user code themselves. Delta only
    // compares fixed components; cyclic only unwraps a sole sym projector.
    // Their descendants remain checked for user callbacks.
    pub(crate) fn intrinsic_normalization_head(symbol: Symbol) -> bool {
        symbol.get_evaluation_info().is_none()
            && (symbol == ETS.metric
                || symbol == ETS.delta
                || symbol == SPENSO_TAG.dot
                || symbol.get_normalization_function().is_none()
                || symbol == AIND_SYMBOLS.dind
                || symbol == *shadowing::CYCLIC)
    }

    /// Index substitution can retain a planned interface when reconstruction
    /// invokes only intrinsic normalization, including within scalar metadata.
    pub fn normalization_is_intrinsic(value: AtomView<'_>) -> bool {
        Self::normalization_is_intrinsic_with_metric(value, &mut None)
    }

    fn intrinsic_head_with_metric(symbol: Symbol, metric: &mut Option<Symbol>) -> bool {
        #[cfg(test)]
        {
            #[cfg(debug_assertions)]
            let debug = DebugInference::is_active();
            #[cfg(not(debug_assertions))]
            let debug = false;
            if debug {
                tests::DEBUG_INTRINSIC_HEAD_CHECKS.with(|count| count.set(count.get() + 1));
            } else {
                tests::INTRINSIC_HEAD_CHECKS.with(|count| count.set(count.get() + 1));
            }
        }
        symbol.get_evaluation_info().is_none()
            && (symbol == *metric.get_or_insert_with(|| ETS.metric)
                || symbol == ETS.delta
                || symbol == SPENSO_TAG.dot
                || symbol.get_normalization_function().is_none()
                || symbol == AIND_SYMBOLS.dind
                || symbol == *shadowing::CYCLIC)
    }

    fn normalization_is_intrinsic_with_metric(
        value: AtomView<'_>,
        metric: &mut Option<Symbol>,
    ) -> bool {
        let mut intrinsic = true;
        let mut slots = SlotMatcher::default();
        value.visitor(&mut |node| {
            if let AtomView::Fun(function) = node {
                intrinsic &= Self::intrinsic_head_with_metric(function.get_symbol(), metric);
                if let SlotMatch::Explicit(slot) = slots.classify(node) {
                    intrinsic &= Self::intrinsic_slot_metadata(slot, metric);
                    return false;
                }
            }
            intrinsic
        });
        intrinsic
    }

    fn intrinsic_slot_metadata(slot: SlotView<'_>, metric: &mut Option<Symbol>) -> bool {
        // Substitution copies a complete index identity. Its payload is not
        // rebuilt, so callbacks written inside it cannot change tensor rank.
        (slot.representation().wrapper().is_none()
            || Self::intrinsic_head_with_metric(slot.representation().head(), metric))
            && Self::normalization_is_intrinsic_with_metric(slot.dimension(), metric)
    }

    fn normalization_is_intrinsic_cached(&mut self, value: AtomView<'_>) -> bool {
        let mut metric = None;
        let mut intrinsic = true;
        value.visitor(&mut |node| {
            if let AtomView::Fun(function) = node {
                if intrinsic
                    && !node.needs_normalization()
                    && node.get_byte_size() <= Self::CACHE_KEY_BYTES
                {
                    if !self.intrinsic_functions.contains(node.get_data()) {
                        intrinsic = Self::normalization_is_intrinsic_with_metric(node, &mut metric);
                        if intrinsic && self.intrinsic_functions.len() < Self::CACHE_ENTRIES {
                            self.intrinsic_functions.insert(node.get_data().to_vec());
                        }
                    }
                    // This complete subtree has already been inspected. Only
                    // positive normalized proofs survive across calls.
                    return false;
                }
                // Preserve the visitor's original failure-tail behavior: false
                // prunes this node, but remaining siblings are still inspected.
                intrinsic &= Self::intrinsic_head_with_metric(function.get_symbol(), &mut metric);
                if let SlotMatch::Explicit(slot) = self.slots.classify(node) {
                    intrinsic &= Self::intrinsic_slot_metadata(slot, &mut metric);
                    return false;
                }
            }
            intrinsic
        });
        intrinsic
    }

    /// Certify that index rewriting cannot invoke a user callback through tensor
    /// leaves or scalar metadata. Only intrinsic metric, delta, dual and cyclic
    /// normalizers are admitted; their descendants remain checked.
    pub fn rewrites_preserve_leaf_interfaces(&mut self, value: AtomView<'_>) -> bool {
        self.normalization_is_intrinsic_cached(value)
            && self.algebra_preserves_leaf_interfaces_impl(value, true)
    }

    /// Scalar algebra treats tensor leaves as indeterminates. Reuse their
    /// interfaces only when inference cannot run callbacks or change power semantics.
    pub fn algebra_preserves_leaf_interfaces(&mut self, value: AtomView<'_>) -> bool {
        self.algebra_preserves_leaf_interfaces_impl(value, false)
    }

    fn algebra_preserves_leaf_interfaces_impl(
        &mut self,
        value: AtomView<'_>,
        intrinsic_proven: bool,
    ) -> bool {
        if matches!(value, AtomView::Fun(_))
            && !value.needs_normalization()
            && value.get_byte_size() <= Self::CACHE_KEY_BYTES
            && self.proven_leaf_interfaces.contains(value.get_data())
        {
            return true;
        }
        match value {
            AtomView::Add(sum) => sum
                .iter()
                .all(|term| self.algebra_preserves_leaf_interfaces_impl(term, intrinsic_proven)),
            AtomView::Mul(product) => product.iter().all(|factor| {
                self.algebra_preserves_leaf_interfaces_impl(factor, intrinsic_proven)
            }),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                // Existing occurrence validation does not multiply power counts.
                // Only opaque dots of compact vectors expose no repeated indices.
                let scalar_dot = matches!(base, AtomView::Fun(function)
                if self.direct_dot_ports(function).is_some_and(|ports| {
                    ports.iter().all(|slot| matches!(slot.aind, PartialIndex::Open(_)))
                }));
                !Self::has_structured_syntax(exponent)
                    && (!Self::has_structured_syntax(base) || scalar_dot)
                    && self.algebra_preserves_leaf_interfaces_impl(base, intrinsic_proven)
                    && self.algebra_preserves_leaf_interfaces_impl(exponent, intrinsic_proven)
            }
            AtomView::Fun(function) if Self::has_structured_syntax(value) => {
                let preserved = self.function_preserves_leaf_interface(function, intrinsic_proven);
                if preserved
                    && !value.needs_normalization()
                    && value.get_byte_size() <= Self::CACHE_KEY_BYTES
                    && self.proven_leaf_interfaces.len() < Self::CACHE_ENTRIES
                {
                    self.proven_leaf_interfaces
                        .insert(value.get_data().to_vec());
                }
                preserved
            }
            // Scalar functions own their metadata, so expansion leaves it opaque.
            _ => true,
        }
    }

    fn function_preserves_leaf_interface(
        &mut self,
        function: symbolica::atom::representation::FunView<'_>,
        intrinsic_proven: bool,
    ) -> bool {
        #[cfg(test)]
        {
            #[cfg(debug_assertions)]
            let debug = DebugInference::is_active();
            #[cfg(not(debug_assertions))]
            let debug = false;
            if debug {
                tests::DEBUG_LEAF_PROOF_CALLS.with(|count| count.set(count.get() + 1));
            } else {
                tests::LEAF_PROOF_CALLS.with(|count| count.set(count.get() + 1));
            }
        }
        let value = function.as_view();
        let symbol = function.get_symbol();
        let arguments = function.iter().collect::<Vec<_>>();
        if self.direct_dot_ports(function).is_some() {
            // Compact dots consume their operands' ports inside the
            // opaque function. Full inference still checks compatibility.
            return true;
        }
        if arguments.iter().all(|argument| {
            self.slots
                .parse::<LibraryRep, AbstractIndex>(*argument)
                .is_ok_and(|slot| {
                    !matches!(slot.aind(), AbstractIndex::Open { .. })
                        && slot.rep().rep.is_self_dual()
                })
        }) && self
            .builtin_tensor_structure(symbol, &arguments, false)
            .is_ok_and(|structure| structure.is_some())
        {
            return true;
        }
        let bound_ports = arguments
            .iter()
            .any(|argument| self.compact_vector_port(*argument));
        let interface = self.direct_leaf_interface(function, self.leaf_inference);
        let interface = if bound_ports || interface.is_none() {
            let plain = |head: Symbol| {
                head.get_normalization_function().is_none()
                    && head.get_evaluation_info().is_none()
                    && !head.is_symmetric()
                    && !head.is_antisymmetric()
                    && !head.is_cyclesymmetric()
                    && !head.is_linear()
                    && !head.is_flat()
            };
            // This proof is also used by scalar algebra without an outer
            // callback guard. Inspect opaque metadata before validation,
            // which must never execute a callback merely to reject reuse.
            if !symbol.has_tag(&SPENSO_TAG.tensor)
                || !plain(symbol)
                || (!intrinsic_proven && !Self::normalization_is_intrinsic(value))
            {
                return false;
            }
            for argument in &arguments {
                if self.compact_vector_port(*argument) {
                    let AtomView::Fun(vector) = argument else {
                        unreachable!()
                    };
                    if !plain(vector.get_symbol())
                        || !self
                            .slots
                            .vector_argument(*vector)
                            .and_then(|port| {
                                self.slots.parse_representation::<LibraryRep>(port).ok()
                            })
                            .is_some_and(|rep| rep.rep.is_self_dual())
                    {
                        return false;
                    }
                }
            }
            if self.validate_tensor_metadata(function).is_err() {
                return false;
            }
            // Scalar algebra preserves these normalized leaves. Read
            // their ports without the constructor's dummy materialization
            // of function-valued metadata, including opaque scalar wrappers.
            self.direct_leaf_interface(function, LeafInference::Observe)
        } else {
            interface
        };
        interface.is_some_and(|interface| {
            interface.logical_slots().iter().all(|slot| {
                matches!(slot.aind, PartialIndex::Explicit(_)) && slot.rep().rep.is_self_dual()
            })
        })
    }

    fn direct_dot_ports(
        &mut self,
        function: symbolica::atom::representation::FunView<'_>,
    ) -> Option<[PartialSlot; 2]> {
        if !(function.get_symbol() == ETS.metric || function.get_symbol() == SPENSO_TAG.dot)
            || function.get_nargs() != 2
        {
            return None;
        }
        let mut ports = Vec::with_capacity(2);
        for operand in function.iter() {
            let AtomView::Fun(operand) = operand else {
                return None;
            };
            if operand
                .iter()
                .any(|argument| self.compact_vector_port(argument))
            {
                return None;
            }
            let interface = self.direct_leaf_interface(operand, self.leaf_inference)?;
            let slots = interface.logical_slots();
            let [slot] = slots.as_slice() else {
                return None;
            };
            ports.push(*slot);
        }
        ports.try_into().ok()
    }

    fn is_builtin_tensor(symbol: Symbol) -> bool {
        symbol == *EPSILON_SYMBOL
            || symbol == ETS.metric
            || symbol == ETS.flat
            || symbol == AGS.gamma
            || symbol == AGS.gamma5
            || symbol == AGS.projm
            || symbol == AGS.projp
            || symbol == AGS.sigma
            || symbol == CS.f
            || symbol == CS.t
    }

    fn builtin_tensor_structure(
        &mut self,
        symbol: Symbol,
        arguments: &[AtomView<'_>],
        logical_arguments: bool,
    ) -> InferenceResult<Option<Canonicalized<ExplicitKey<AbstractIndex>>>> {
        if !Self::is_builtin_tensor(symbol) {
            return Ok(None);
        }

        let compact_metric = symbol == ETS.metric
            && arguments.len() == 2
            && arguments.iter().all(|argument| {
                self.slots
                    .parse::<LibraryRep, AbstractIndex>(*argument)
                    .is_err()
                    && self
                        .slots
                        .parse_representation::<LibraryRep>(*argument)
                        .is_err()
                    && !Self::is_chain_placeholder(*argument)
                    && Self::has_structured_syntax(*argument)
            });
        if compact_metric {
            return Ok(None);
        }

        let ports = arguments
            .iter()
            .map(|argument| {
                if Self::is_chain_placeholder(*argument) {
                    Ok(None)
                } else {
                    self.slots
                        .parse::<LibraryRep, AbstractIndex>(*argument)
                        .map(|slot| Some(slot.rep()))
                        .or_else(|_| {
                            self.slots
                                .parse_representation::<LibraryRep>(*argument)
                                .map(Some)
                        })
                        .map_err(|_| ())
                        .or_else(|_| {
                            // A consumed compact vector contributes its declared port,
                            // not the result of constructing a fresh component. Replaying
                            // its normalizer here can change rank before the rewrite runs.
                            if let AtomView::Fun(vector) = *argument
                                && self.compact_vector_port(*argument)
                            {
                                self.validate_tensor_metadata(vector).map_err(|_| ())?;
                                let compact = self.slots.vector_argument(vector).ok_or(())?;
                                return self
                                    .slots
                                    .parse_representation::<LibraryRep>(compact)
                                    .map(Some)
                                    .map_err(|_| ());
                            }
                            // Consumed operands have a logical compact axis even
                            // when the enclosing query observes graph storage.
                            // Keep their logical interface cache separate; only
                            // the mode-independent slot matcher is shared.
                            let mut operand = Self {
                                slots: mem::take(&mut self.slots),
                                leaf_inference: LeafInference::Observe,
                                ..Self::default()
                            };
                            let interface = operand.infer_validated(*argument);
                            self.slots = operand.slots;
                            let interface = interface.map_err(|_| ())?;
                            let slots = interface.logical_slots();
                            match slots.as_slice() {
                                [slot] if matches!(slot.aind, PartialIndex::Open(_)) => {
                                    Ok(Some(slot.rep()))
                                }
                                _ => Err(()),
                            }
                        })
                }
            })
            .collect::<Result<Vec<_>, _>>();
        let dimension = |representation: LibraryRep| {
            ports
                .as_ref()
                .ok()
                .and_then(|ports| {
                    ports
                        .iter()
                        .flatten()
                        .find(|port| port.rep == representation)
                })
                .map(|port| port.dim)
        };
        let invalid = |factory: &'static str, signature: &'static str| {
            TensorInferenceError::InvalidBuiltinSignature { factory, signature }
        };

        let (factory, signature, expected): (
            _,
            _,
            Option<Canonicalized<ExplicitKey<AbstractIndex>>>,
        ) = if symbol == ETS.metric || symbol == ETS.flat {
            let factory = if symbol == ETS.metric { "g" } else { "flat" };
            let ports = ports
                .as_ref()
                .map_err(|_| invalid(factory, "two compatible representation ports"))?;
            if arguments.len() != 2 {
                return Err(invalid(factory, "two compatible representation ports"));
            }
            let visible = ports.iter().flatten().copied().collect::<Vec<_>>();
            let compatible = visible.len() < 2
                || if symbol == ETS.metric {
                    visible[0] == visible[1] || visible[0].matches(&visible[1])
                } else {
                    visible[0] == visible[1]
                };
            if !compatible {
                return Err(invalid(factory, "two compatible representation ports"));
            }
            let expected = match visible.as_slice() {
                [left, right] => Some(ExplicitKey::from_iter([*left, *right], symbol, None)),
                _ => None,
            };
            (factory, "two compatible representation ports", expected)
        } else if symbol == *EPSILON_SYMBOL {
            let signature = "equal representation ports or contracted vectors";
            let ports = ports.as_ref().map_err(|_| invalid("epsilon", signature))?;
            let visible = ports
                .iter()
                .copied()
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| invalid("epsilon", signature))?;
            if visible.is_empty() || visible.iter().any(|port| *port != visible[0]) {
                return Err(invalid("epsilon", signature));
            }
            (
                "epsilon",
                signature,
                Some(ExplicitKey::from_iter(visible, symbol, None)),
            )
        } else if symbol == AGS.gamma {
            (
                "gamma",
                if logical_arguments {
                    "one Minkowski and two equal bispinor ports"
                } else {
                    "one Minkowski and two four-dimensional bispinor ports"
                },
                Some(
                    AGS.gamma_strct(
                        dimension(Minkowski {}.into()).unwrap_or(Dimension::Concrete(1)),
                    ),
                ),
            )
        } else if symbol == AGS.gamma5 || symbol == AGS.projm || symbol == AGS.projp {
            let dimension = dimension(Bispinor {}.into()).unwrap_or(Dimension::Concrete(1));
            let (factory, expected) = if symbol == AGS.gamma5 {
                ("gamma5", AGS.gamma5_strct(dimension))
            } else if symbol == AGS.projm {
                ("projm", AGS.projm_strct(dimension))
            } else {
                ("projp", AGS.projp_strct(dimension))
            };
            (factory, "two equal bispinor ports", Some(expected))
        } else if symbol == AGS.sigma {
            (
                "sigma",
                if logical_arguments {
                    "two equal Minkowski and two equal bispinor ports"
                } else {
                    "two equal Minkowski and two four-dimensional bispinor ports"
                },
                Some(
                    AGS.sigma_strct(
                        dimension(Minkowski {}.into()).unwrap_or(Dimension::Concrete(1)),
                    ),
                ),
            )
        } else if symbol == CS.f {
            (
                "f",
                "three equal adjoint ports",
                Some(
                    CS.f_strct(dimension(ColorAdjoint {}.into()).unwrap_or(Dimension::Concrete(1))),
                ),
            )
        } else if symbol == CS.t {
            let adjoint_dimension =
                dimension(ColorAdjoint {}.into()).unwrap_or(Dimension::Concrete(1));
            let fundamental_dimension = dimension(ColorFundamental {}.into())
                .or_else(|| dimension(ColorAntiFundamental {}.into()))
                .unwrap_or(Dimension::Concrete(1));
            (
                "t",
                "adjoint, fundamental, and antifundamental ports with matching color dimensions",
                Some(CS.t_strct(fundamental_dimension, adjoint_dimension)),
            )
        } else {
            unreachable!("every predefined tensor symbol was handled")
        };

        let ports = ports.map_err(|_| invalid(factory, signature))?;
        let placeholders = arguments
            .iter()
            .enumerate()
            .filter_map(|(position, argument)| {
                let AtomView::Var(variable) = *argument else {
                    return None;
                };
                let symbol = variable.get_symbol();
                Self::is_chain_placeholder(*argument).then_some((position, symbol))
            })
            .collect::<Vec<_>>();
        if !placeholders.is_empty()
            && (!matches!(placeholders.as_slice(), [(_, input), (_, output)] if *input != *output)
                || !placeholders
                    .iter()
                    .any(|(_, symbol)| *symbol == SPENSO_TAG.chain_in)
                || !placeholders
                    .iter()
                    .any(|(_, symbol)| *symbol == SPENSO_TAG.chain_out))
        {
            return Err(invalid(factory, signature));
        }
        let Some(expected) = expected else {
            return Ok(None);
        };
        let expected_ports = expected
            .canonical()
            .external_reps_iter()
            .collect::<Vec<_>>();
        // A metric signature is built from its observed ports; its symmetric
        // Atom argument order need not match representation storage order.
        // Standalone factory atoms still spell canonical storage ports, while
        // chain/trace words retain the factory's logical order around in/out.
        let mut expected_ports =
            if symbol == ETS.metric || logical_arguments || !placeholders.is_empty() {
                expected.layout().canonical_to_logical(&expected_ports)
            } else {
                expected_ports
            };
        // A bound formal word inherits its spin space from the chain/trace
        // channel. The factory's bis(4) default remains the standalone tensor
        // signature; four-dimensional algebra identities have their own guards.
        if logical_arguments
            && (symbol == AGS.gamma || symbol == AGS.sigma)
            && let Some(spin) = ports
                .iter()
                .flatten()
                .find(|port| port.rep == LibraryRep::from(Bispinor {}))
        {
            for expected in &mut expected_ports {
                if expected.rep == spin.rep {
                    expected.dim = spin.dim;
                }
            }
        }
        if ports.len() != expected_ports.len()
            || ports
                .iter()
                .zip(&expected_ports)
                .any(|(actual, expected)| actual.is_some_and(|actual| actual != *expected))
        {
            return Err(invalid(factory, signature));
        }
        if !placeholders.is_empty() {
            let input = placeholders
                .iter()
                .find(|(_, symbol)| *symbol == SPENSO_TAG.chain_in)
                .expect("placeholder names were validated")
                .0;
            let output = placeholders
                .iter()
                .find(|(_, symbol)| *symbol == SPENSO_TAG.chain_out)
                .expect("placeholder names were validated")
                .0;
            let input = expected_ports[input];
            let output = expected_ports[output];
            let oriented =
                input.rep.is_self_dual() || (input.rep.is_base() && output.rep.is_dual());
            if !input.matches(&output) || !oriented {
                return Err(invalid(factory, signature));
            }
        }
        Ok(Some(expected))
    }

    fn validate_builtin_placeholder_channels(
        &mut self,
        value: AtomView<'_>,
        input: Representation<LibraryRep>,
        output: Representation<LibraryRep>,
    ) -> InferenceResult<()> {
        match value {
            AtomView::Add(add) => {
                for term in add.iter() {
                    self.validate_builtin_placeholder_channels(term, input, output)?;
                }
            }
            AtomView::Mul(mul) => {
                for factor in mul.iter() {
                    self.validate_builtin_placeholder_channels(factor, input, output)?;
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                self.validate_builtin_placeholder_channels(base, input, output)?;
                self.validate_builtin_placeholder_channels(exponent, input, output)?;
            }
            AtomView::Fun(function) => {
                let arguments = function.iter().collect::<Vec<_>>();
                if arguments
                    .iter()
                    .any(|argument| Self::is_chain_placeholder(*argument))
                {
                    let resolved = arguments
                        .iter()
                        .map(|argument| match *argument {
                            AtomView::Var(variable)
                                if variable.get_symbol() == SPENSO_TAG.chain_in =>
                            {
                                input.to_symbolic([])
                            }
                            AtomView::Var(variable)
                                if variable.get_symbol() == SPENSO_TAG.chain_out =>
                            {
                                output.to_symbolic([])
                            }
                            argument => argument.to_owned(),
                        })
                        .collect::<Vec<_>>();
                    let resolved = resolved.iter().map(Atom::as_view).collect::<Vec<_>>();
                    self.builtin_tensor_structure(function.get_symbol(), &resolved, true)?;
                }
                for argument in arguments {
                    if !Self::is_chain_placeholder(argument) {
                        self.validate_builtin_placeholder_channels(argument, input, output)?;
                    }
                }
            }
            _ => {}
        }
        Ok(())
    }
}

impl InterfaceInference {
    fn publish_admission(
        &mut self,
        expression: Atom,
        structure: PartialStructure,
        occurrences: ExplicitIndexOccurrences,
    ) -> SymbolicTensor<PartialStructure> {
        self.purpose = ObservationPurpose::Validate;
        let structure = SymbolicTensor::canonical_structure(structure);
        let leaves = self.admitted_interfaces.declare_root(
            expression.as_view(),
            structure.clone(),
            Arc::new(Mutex::new(mem::take(&mut self.admitted_leaf_interfaces))),
        );
        SymbolicTensor::from_published_parts(
            expression,
            structure,
            Arc::new(occurrences),
            OnceLock::from(Arc::new(mem::take(&mut self.admitted_interfaces))),
            leaves,
        )
    }

    pub fn infer_validated(&mut self, atom: AtomView<'_>) -> InferenceResult<PartialStructure> {
        self.infer_observed(atom).map(|(interface, _)| interface)
    }

    fn infer_observed(
        &mut self,
        atom: AtomView<'_>,
    ) -> InferenceResult<(PartialStructure, ExplicitIndexOccurrences)> {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::InterfaceDiscovery,
        );
        #[cfg(test)]
        {
            #[cfg(debug_assertions)]
            let debug = DebugInference::is_active();
            #[cfg(not(debug_assertions))]
            let debug = false;
            if debug {
                tests::DEBUG_INFERENCE_CALLS.with(|count| count.set(count.get() + 1));
            } else {
                tests::INFERENCE_CALLS.with(|count| count.set(count.get() + 1));
            }
        }
        if self.purpose == ObservationPurpose::Publish {
            self.admitted_interfaces = AdmittedInterfaces::for_expression(atom);
        }
        let occurrences = Self::scoped_occurrences(atom, &mut self.slots)?;
        self.validate_tensor_syntax_after_scopes(atom)?;
        occurrences.validate()?;
        if self.leaf_inference == LeafInference::Materialize {
            self.dummies.reserve_indices(atom);
        }
        // Scalar expressions still pass the same scope and signature checks;
        // Scalar-tagged metadata does not expose tensor ports.
        let interface = if Self::has_structured_syntax(atom) {
            self.infer_view(atom)?
        } else {
            PartialStructure::from_logical_slots([])
        };
        Ok((interface, occurrences))
    }

    fn validate_tensor_syntax_after_scopes(&mut self, atom: AtomView<'_>) -> InferenceResult<()> {
        let mut syntax_error = None;
        let mut checked_functions = AHashSet::new();
        atom.visitor(&mut |value| {
        if syntax_error.is_some() {
            return false;
        }

        match value {
            AtomView::Var(variable) if variable.get_symbol().has_tag(&SPENSO_TAG.tensor) => {
                syntax_error = Some(TensorInferenceError::invalid(format!(
                    "tensor symbol `{}` must be called with at least one structural port",
                    variable.get_symbol()
                )));
            }
            AtomView::Fun(function) => {
                if value.get_byte_size() <= InterfaceInference::CACHE_KEY_BYTES
                    && checked_functions.contains(&value)
                {
                    return false;
                }
                let symbol = function.get_symbol();
                // Explicit scalar functions own their metadata, which may itself
                // contain tensor syntax without exposing any tensor ports.
                if SPENSO_TAG.is_scalar_metadata(symbol)
                    || symbol == SPENSO_TAG.scalar
                    || symbol == SPENSO_TAG.pure_scalar
                {
                    return false;
                }
                let arguments = function.iter().collect::<Vec<_>>();
                // A representation carrying an index must be a valid slot. Falling back
                // to its bare representation would turn contracted indices into open ports.
                if symbol.has_tag(&SPENSO_TAG.representation)
                    && arguments.len() > 1
                    && (arguments.len() != 2
                        || self.slots.parse::<LibraryRep, AbstractIndex>(value).is_err())
                {
                    syntax_error = Some(TensorInferenceError::invalid(format!(
                        "invalid indexed representation `{value}`; cook nested indices before inferring a tensor interface"
                    )));
                    return false;
                }
                // A complete slot/representation owns its dimension and index
                // metadata. A valid symbol index can also name a tensor head;
                // it is not a standalone tensor expression at this boundary.
                if self.slots.parse::<LibraryRep, AbstractIndex>(value).is_ok()
                    || self.slots.parse_representation::<LibraryRep>(value).is_ok()
                {
                    return false;
                }
                if symbol == SPENSO_TAG.bracket && arguments.is_empty() {
                    syntax_error = Some(TensorInferenceError::invalid("ordered tensor products require at least one operand"));
                    return false;
                }
                if symbol == SPENSO_TAG.dot && arguments.len() != 2 {
                    syntax_error = Some(TensorInferenceError::invalid("dot requires exactly two operands"));
                    return false;
                }
                if symbol == SPENSO_TAG.chain && arguments.len() < 2 {
                    syntax_error = Some(TensorInferenceError::invalid("chain requires explicit start and end ports"));
                    return false;
                }
                if symbol == SPENSO_TAG.trace {
                    let Some(representation) = arguments.first() else {
                        syntax_error = Some(TensorInferenceError::invalid("trace requires a representation"));
                        return false;
                    };
                    if self.slots.parse_representation::<LibraryRep>(*representation).is_err() {
                        syntax_error =
                            Some(TensorInferenceError::invalid("trace metadata is not a Spenso representation"));
                        return false;
                    }
                }
                if symbol.has_tag(&SPENSO_TAG.broadcast) && arguments.len() != 1 {
                    syntax_error = Some(TensorInferenceError::invalid(format!(
                        "broadcast function `{symbol}` requires exactly one argument"
                    )));
                    return false;
                }
                if !symbol.has_tag(&SPENSO_TAG.tensor) {
                    // Typed admission uses the same tagged tensor grammar as
                    // its graph consumers. Ordinary scalar functions remain
                    // opaque, but cannot silently discard direct tensor ports.
                    if !composition::is_composite_head(symbol, &SPENSO_TAG)
                        && symbol != SPENSO_TAG.dot
                        && arguments.iter().any(|argument| {
                            self.slots.parse::<LibraryRep, AbstractIndex>(*argument).is_ok()
                                || self.slots.parse_representation::<LibraryRep>(*argument).is_ok()
                        })
                    {
                        syntax_error = Some(TensorInferenceError::invalid(format!(
                            "function `{symbol}` carries structural ports but is not tagged as a tensor"
                        )));
                        return false;
                    }
                    if self.purpose == ObservationPurpose::Publish
                        && !Self::has_structured_syntax(value)
                        && !composition::is_composite_head(symbol, &SPENSO_TAG)
                        && symbol != SPENSO_TAG.dot
                    {
                        // A product can flatten or acquire a sign and expose this
                        // scalar leaf independently of its originally admitted root.
                        self.admitted_interfaces.insert(value, PartialStructure::from_logical_slots([]));
                        self.admitted_leaf_interfaces.insert(value.to_owned(), (false, Vec::new()));
                    }
                    return true;
                }

                let compact_metric = symbol == ETS.metric
                    && arguments.len() == 2
                    && arguments.iter().all(|argument| {
                        self.slots.parse::<LibraryRep, AbstractIndex>(*argument).is_err()
                            && self.slots.parse_representation::<LibraryRep>(*argument).is_err()
                            && Self::has_structured_syntax(*argument)
                    });
                let builtin = match self.builtin_tensor_structure(symbol, &arguments, false) {
                    Ok(builtin) => builtin,
                    Err(error) => {
                        syntax_error = Some(error);
                        return false;
                    }
                };

                if let Some(structure) = &builtin
                    && !symbol.has_tag(&SPENSO_TAG.rank1)
                    && arguments.iter().all(|argument| self.slots.parse::<LibraryRep, AbstractIndex>(*argument)
                        .is_ok_and(|slot| !matches!(slot.aind(), AbstractIndex::Open { .. })))
                {
                    // The builtin signature established every explicit port.
                    // Slots own their metadata, so neither generic leaf ordering
                    // nor another child syntax walk can add a validation here.
                    if checked_functions.len() < InterfaceInference::CACHE_ENTRIES
                        && value.get_byte_size() <= InterfaceInference::CACHE_KEY_BYTES
                    {
                        self.builtin_interface(value, &arguments, structure);
                        checked_functions.insert(value);
                    }
                    return false;
                }

                let ports = arguments
                    .iter()
                    .enumerate()
                    .filter(|(position, argument)| {
                        if symbol.has_tag(&SPENSO_TAG.rank1)
                            && let Some(component) = self.slots.concrete_component(**argument)
                        {
                            if component.is_err() {
                                syntax_error = Some(TensorInferenceError::invalid(format!(
                                    "rank-one tensor function `{symbol}` requires one nonnegative integer concrete component"
                                )));
                            }
                            // δ(cind(k), port) is a basis vector. Its first
                            // concrete marker selects the component; only the
                            // final argument is a tensor port.
                            if symbol == ETS.delta
                                && *position == 0
                                && arguments.len() == 2
                                && matches!(**argument, AtomView::Fun(marker)
                                    if marker.get_symbol() == AIND_SYMBOLS.cind)
                            {
                                return false;
                            }
                            // This counts the declared rank for syntax only. The
                            // existing leaf reader exposes no axis for a fixed component.
                            return true;
                        }
                        self.slots.parse::<LibraryRep, AbstractIndex>(**argument).is_ok()
                            || self.slots.parse_representation::<LibraryRep>(**argument).is_ok()
                            || (!symbol.has_tag(&SPENSO_TAG.rank1)
                                && self.compact_vector_port(**argument))
                            || matches!(
                                **argument,
                                AtomView::Var(variable)
                                    if variable.get_symbol() == SPENSO_TAG.chain_in
                                        || variable.get_symbol() == SPENSO_TAG.chain_out
                            )
                    })
                    .map(|(position, _)| position)
                    .collect::<Vec<_>>();
                if syntax_error.is_some() {
                    return false;
                }
                let ports_are_final = ports
                    .iter()
                    .copied()
                    .eq(arguments.len().saturating_sub(ports.len())..arguments.len());
                if !ports_are_final && !compact_metric && builtin.is_none() {
                    syntax_error = Some(TensorInferenceError::invalid(format!(
                        "tensor function `{symbol}` requires scalar arguments before structural ports"
                    )));
                } else if symbol.has_tag(&SPENSO_TAG.rank1) && ports.len() != 1 {
                    syntax_error = Some(TensorInferenceError::invalid(format!(
                        "rank-one tensor function `{symbol}` requires exactly one final structural port"
                    )));
                }
                // Explicit built-ins and directly readable leaves have no
                // materialization effects during syntax validation. The first
                // visit still checks every child before a later occurrence is skipped.
                if syntax_error.is_none()
                    && checked_functions.len() < InterfaceInference::CACHE_ENTRIES
                    && value.get_byte_size() <= InterfaceInference::CACHE_KEY_BYTES
                    && self.direct_leaf_interface(function, self.leaf_inference).is_some()
                {
                    checked_functions.insert(value);
                }
            }
            _ => {}
        }
        true
    });
        if let Some(error) = syntax_error {
            return Err(error);
        }

        Ok(())
    }
}

#[derive(Default)]
pub struct InterfaceInference {
    purpose: ObservationPurpose,
    admitted_leaf_interfaces: HashMap<Atom, (bool, Vec<Atom>)>,
    admitted_interfaces: AdmittedInterfaces,
    reusable_interfaces: AHashMap<Vec<u8>, PartialStructure>,
    // Shape inference alone does not establish this stronger rewrite proof.
    proven_leaf_interfaces: AHashSet<Vec<u8>>,
    // A leaf-interface proof can leave scalar metadata opaque; this separate
    // certificate covers every function in a complete normalized subtree.
    intrinsic_functions: AHashSet<Vec<u8>>,
    slots: SlotMatcher,
    dummies: ParseState<AbstractIndex>,
    leaf_inference: LeafInference,
}

/// Validation and admission run the same syntax and multiplicity checks.
/// Only admission keeps a carrier's regional interfaces and written leaf ports.
#[derive(Clone, Copy, Default, PartialEq, Eq)]
enum ObservationPurpose {
    #[default]
    Validate,
    Publish,
}

#[derive(Clone, Copy, Default, PartialEq)]
pub(super) enum LeafInference {
    #[default]
    Materialize,
    Observe,
    ObserveEncoded,
    // A parsed leaf keeps locally sewn axes in its graph storage.
    ObserveStorage,
}

impl LeafInference {
    fn for_declaration(interface: &PartialStructure) -> Self {
        if interface.logical_slots().iter().any(|slot| {
            matches!(
                slot.aind,
                PartialIndex::Explicit(AbstractIndex::Open { .. })
            )
        }) {
            Self::ObserveEncoded
        } else {
            Self::Observe
        }
    }
}

impl InterfaceInference {
    /// Observe local replacement ports without materialization or loss of
    /// encoded open-index owners. Scalar RHS expressions also undergo the
    /// normal syntax checks, including rejection of bare tagged tensor names.
    ///
    /// Observation does not rebuild the expression or invoke its normalization
    /// callbacks. Native scalar metadata remains opaque, and encoded open
    /// owners and axis identities remain explicit in the returned interface.
    /// A closed interface alone does not certify ordinary commuting algebra.
    /// Since this observes an expression rather than an established carrier,
    /// an exact zero has a scalar interface; retain a typed zero's original
    /// interface separately when checking a transformation.
    ///
    /// # Errors
    /// Rejects invalid tensor syntax, malformed structural ports, incompatible
    /// summand interfaces, and excessive explicit-index occurrences.
    ///
    /// # Example
    ///
    /// ```
    /// use idenso::tensor::inference::{InterfaceInference, TensorInferenceError};
    /// use spenso::structure::partial::PartialStructureExt;
    /// use symbolica::atom::{Atom, AtomCore};
    /// let expression = Atom::num(7);
    /// let interface = InterfaceInference::replacement_interface(expression.as_view())?;
    /// assert!(interface.logical_slots().is_empty());
    /// # Ok::<(), TensorInferenceError>(())
    /// ```
    pub fn replacement_interface(atom: AtomView<'_>) -> InferenceResult<PartialStructure> {
        Self::replacement_observation(atom, &mut SlotMatcher::default())
            .map(|(interface, _)| interface)
    }

    pub(super) fn replacement_observation(
        atom: AtomView<'_>,
        slots: &mut SlotMatcher,
    ) -> InferenceResult<(PartialStructure, ExplicitIndexOccurrences)> {
        let occurrences = Self::scoped_occurrences(atom, slots)?;
        let mut inference = Self {
            leaf_inference: LeafInference::ObserveEncoded,
            ..Self::default()
        };
        inference.validate_tensor_syntax_after_scopes(atom)?;
        let interface = if Self::has_structured_syntax(atom) {
            Self::merge_explicit_interface_sequence(&[inference.infer_view(atom)?])?
                .canonicalize_open_ports()
        } else {
            PartialStructure::from_logical_slots([])
        };
        Ok((interface, occurrences))
    }

    const CACHE_ENTRIES: usize = 256;
    const CACHE_KEY_BYTES: usize = 256;

    pub(super) fn partial_index(
        index: AbstractIndex,
        policy: LeafInference,
    ) -> PartialIndex<AbstractIndex> {
        match index {
            AbstractIndex::Open { axis, .. }
                if !matches!(
                    policy,
                    LeafInference::ObserveEncoded | LeafInference::ObserveStorage
                ) =>
            {
                PartialIndex::open(axis)
            }
            index => PartialIndex::Explicit(index),
        }
    }

    /// Call only for interfaces independent of materialization callbacks and
    /// fresh dummy identities. Open ports remain occurrence-local positions.
    fn cache_interface(&mut self, atom: AtomView<'_>, interface: &PartialStructure) {
        if self.reusable_interfaces.len() < Self::CACHE_ENTRIES
            && atom.get_byte_size() <= Self::CACHE_KEY_BYTES
        {
            self.reusable_interfaces
                .insert(atom.get_data().to_vec(), interface.clone());
        }
    }

    #[cfg(test)]
    fn infer(&mut self, atom: &Atom) -> InferenceResult<PartialStructure> {
        self.dummies.reserve_indices(atom.as_view());
        self.infer_view(atom.as_view())
    }

    /// A direct compact tagged vector consumes one tensor port. Its scalar
    /// metadata is checked separately; it contributes no external index.
    pub fn compact_vector_port(&mut self, value: AtomView<'_>) -> bool {
        let AtomView::Fun(function) = value else {
            return false;
        };
        function.get_symbol().has_tag(&SPENSO_TAG.rank1)
            && self
                .slots
                .vector_argument(function)
                .is_some_and(|argument| {
                    self.slots
                        .parse_representation::<LibraryRep>(argument)
                        .is_ok()
                })
    }

    /// Validate metadata without materializing a consumed compact vector.
    /// Scalar wrappers keep their contents opaque; an unwrapped tensor must
    /// either be scalar or occupy a canonical consumed port of a generic leaf.
    fn validate_tensor_metadata(
        &mut self,
        function: symbolica::atom::representation::FunView<'_>,
    ) -> InferenceResult<()> {
        let symbol = function.get_symbol();
        for argument in function.iter() {
            if self
                .slots
                .parse::<LibraryRep, AbstractIndex>(argument)
                .is_ok()
                || self
                    .slots
                    .parse_representation::<LibraryRep>(argument)
                    .is_ok()
                || Self::is_chain_placeholder(argument)
                || !Self::has_structured_syntax(argument)
            {
                continue;
            }
            if !symbol.has_tag(&SPENSO_TAG.rank1) && self.compact_vector_port(argument) {
                let AtomView::Fun(vector) = argument else {
                    unreachable!()
                };
                // The vector's earlier arguments remain scalar metadata, not
                // further consumed ports. Its rank-one tag enforces that rule.
                self.validate_tensor_metadata(vector)?;
            } else if !self.infer_view(argument)?.canonical().is_scalar() {
                return Err(TensorInferenceError::invalid(format!(
                    "tensor metadata for `{symbol}` must be scalar"
                )));
            }
        }
        Ok(())
    }

    /// Read ordinary leaf ports in their existing syntactic order. Compact
    /// representations denote occurrence-local open ports; no dummy atom or
    /// function reconstruction is needed to discover their representation.
    fn direct_leaf_interface(
        &mut self,
        function: symbolica::atom::representation::FunView<'_>,
        policy: LeafInference,
    ) -> Option<PartialStructure> {
        let symbol = function.get_symbol();
        let bound_ports = symbol.has_tag(&SPENSO_TAG.tensor)
            && !symbol.has_tag(&SPENSO_TAG.rank1)
            && function
                .iter()
                .any(|argument| self.compact_vector_port(argument));
        if !matches!(self.slots.classify(function.as_view()), SlotMatch::Other)
            || (!bound_ports
                && policy == LeafInference::Materialize
                && (symbol.get_normalization_function().is_some()
                    || symbol.get_evaluation_info().is_some()
                    || symbol.is_symmetric()
                    || symbol.is_antisymmetric()
                    || symbol.is_cyclesymmetric()
                    || symbol.is_linear()))
        {
            return None;
        }
        let mut logical = Vec::new();
        let mut written_slots = Vec::new();
        let mut placeholders = 0;
        let mut explicit_ports = 0;
        let mut function_metadata = false;
        for (position, argument) in function.iter().enumerate() {
            // In d_R(R, a, b, ...), the first argument declares the
            // representation whose generators define the invariant.
            if symbol == CS.d && function.get_nargs() >= 4 && position == 0 {
                continue;
            }
            if Self::is_chain_placeholder(argument) {
                placeholders += 1;
            } else if let Ok(slot) = self.slots.parse::<LibraryRep, AbstractIndex>(argument) {
                explicit_ports += 1;
                let index = Self::partial_index(slot.aind(), policy);
                logical.push(slot.rep().slot(index));
            } else if let Ok(rep) = self.slots.parse_representation::<LibraryRep>(argument) {
                logical.push(rep.slot(PartialIndex::open(logical.len())));
            } else if bound_ports {
                // The normalized bound leaf has already passed metadata checks.
                // Read its surviving direct ports without introducing dummy
                // indices inside the consumed vector or replaying callbacks.
                continue;
            } else {
                // Nested function metadata may contain compact representations,
                // bundles, or callbacks affected by legacy materialization.
                // Keep that existing owner for those forms.
                let mut contains_function = false;
                argument.visitor(&mut |node| {
                    contains_function |= matches!(node, AtomView::Fun(_));
                    !contains_function
                });
                if contains_function && policy == LeafInference::Materialize {
                    return None;
                }
                function_metadata |= contains_function;
                continue;
            }
            if self.purpose == ObservationPurpose::Publish {
                written_slots.push(argument.to_owned());
            }
        }
        if function_metadata {
            // Observation may ignore scalar metadata, but must not hide a
            // structural bundle that the existing fast leaf parser exposes.
            let inferred = match function
                .as_view()
                .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
            {
                Ok(inferred) => inferred.canonical().order(),
                Err(StructureError::EmptyStructure(_)) => 0,
                Err(_) => return None,
            };
            if inferred != explicit_ports {
                return None;
            }
        }
        if self.purpose == ObservationPurpose::Publish
            && written_slots.len() == logical.len() + placeholders
        {
            self.admitted_leaf_interfaces
                .insert(function.as_view().to_owned(), (true, written_slots));
        }
        Some(PartialStructure::from_logical_slots(logical))
    }

    /// Reuse a checked builtin signature without constructing its factory again.
    /// Fully explicit ports need no callback or fresh dummy materialization.
    fn builtin_interface(
        &mut self,
        atom: AtomView<'_>,
        arguments: &[AtomView<'_>],
        structure: &Canonicalized<ExplicitKey<AbstractIndex>>,
    ) -> PartialStructure {
        let canonical_ports = arguments
            .iter()
            .enumerate()
            .map(|(position, argument)| {
                if let Ok(slot) = self.slots.parse::<LibraryRep, AbstractIndex>(*argument) {
                    let index = Self::partial_index(slot.aind(), self.leaf_inference);
                    Some(slot.rep().slot(index))
                } else {
                    self.slots
                        .parse_representation::<LibraryRep>(*argument)
                        .ok()
                        .map(|representation| representation.slot(PartialIndex::open(position)))
                }
            })
            .collect::<Vec<_>>();
        // Only fully explicit built-in leaves are independent of fresh
        // port materialization and normalization callbacks.
        let reusable = canonical_ports.iter().all(
            |port| matches!(port, Some(slot) if matches!(slot.aind, PartialIndex::Explicit(_))),
        );
        // Compact arguments are contracted vectors, and chain placeholders
        // are wiring labels. Neither exposes an external tensor port.
        let logical_ports = if arguments
            .iter()
            .any(|argument| Self::is_chain_placeholder(*argument))
        {
            canonical_ports
        } else {
            structure.layout().canonical_to_logical(&canonical_ports)
        };
        let interface = PartialStructure::from_logical_slots(logical_ports.into_iter().flatten());
        if self.purpose == ObservationPurpose::Publish {
            let written = arguments
                .iter()
                .filter_map(|argument| {
                    if Self::is_chain_placeholder(*argument)
                        || self
                            .slots
                            .parse_representation::<LibraryRep>(*argument)
                            .is_ok()
                    {
                        Some(argument.to_owned())
                    } else {
                        self.slots
                            .parse::<LibraryRep, AbstractIndex>(*argument)
                            .ok()
                            .map(|_| argument.to_owned())
                    }
                })
                .collect();
            // The checked signature owns all argument roles. Consumed tensor
            // operands expose no storage axis; direct ports keep their actual
            // normalized argument order, independently of the factory layout.
            self.admitted_leaf_interfaces
                .insert(atom.to_owned(), (true, written));
        }
        if reusable {
            // Keep repeated indices and logical ordering intact: enclosing
            // products still own contraction and multiplicity validation.
            self.cache_interface(atom, &interface);
            if !atom.needs_normalization()
                && atom.get_byte_size() <= Self::CACHE_KEY_BYTES
                && self.proven_leaf_interfaces.len() < Self::CACHE_ENTRIES
                && interface
                    .logical_slots()
                    .iter()
                    .all(|slot| slot.rep().rep.is_self_dual()
                        && matches!(slot.aind, PartialIndex::Explicit(index) if !matches!(index, AbstractIndex::Open { .. })))
            {
                self.proven_leaf_interfaces.insert(atom.get_data().to_vec());
            }
        }
        interface
    }

    fn cached_interface(&self, atom: AtomView<'_>) -> Option<&PartialStructure> {
        if matches!(atom, AtomView::Fun(_)) && atom.get_byte_size() <= Self::CACHE_KEY_BYTES {
            self.reusable_interfaces.get(atom.get_data())
        } else {
            None
        }
    }

    /// Compose a chain's endpoints and factor spectators. A typed producer
    /// supplies factor interfaces, while raw admission observes its factors.
    fn chain_interface_with<'a>(
        &mut self,
        arguments: &[AtomView<'a>],
        mut factor_interface: impl FnMut(&mut Self, AtomView<'a>) -> InferenceResult<PartialStructure>,
    ) -> InferenceResult<PartialStructure> {
        if arguments.len() < 2 {
            return Err(TensorInferenceError::invalid(
                "chain requires explicit start and end ports",
            ));
        }
        let mut exposed = Vec::new();
        let endpoints = arguments[..2]
            .iter()
            .map(|argument| {
                if let Ok(slot) = self.slots.parse::<LibraryRep, AbstractIndex>(*argument) {
                    let slot = slot
                        .rep()
                        .slot(Self::partial_index(slot.aind(), self.leaf_inference));
                    exposed.push(slot);
                    Ok(slot.rep())
                } else if let Ok(representation) =
                    self.slots.parse_representation::<LibraryRep>(*argument)
                {
                    exposed.push(representation.slot(PartialIndex::open(0)));
                    Ok(representation)
                } else if let AtomView::Fun(vector) = *argument
                    && self.compact_vector_port(*argument)
                {
                    self.validate_tensor_metadata(vector)?;
                    let (rep, dimension) = self.slots.port_representation(*argument).unwrap();
                    Ok(Representation {
                        rep,
                        dim: Dimension::try_from(dimension).map_err(|error| {
                            TensorInferenceError::invalid(error.to_string())
                        })?,
                    })
                } else {
                    Err(TensorInferenceError::invalid(
                        "chain endpoints must be Spenso slots, representations, or compact rank-one tensors",
                    ))
                }
            })
            .collect::<InferenceResult<Vec<_>>>()?;
        let input = endpoints[0];
        let output = endpoints[1];
        if !input.matches(&output)
            || !(input.rep.is_self_dual() || (input.rep.is_base() && output.rep.is_dual()))
        {
            return Err(TensorInferenceError::invalid(
                "chain endpoints do not form a compatible input-to-output channel",
            ));
        }
        let mut interfaces = Vec::new();
        for factor in &arguments[2..] {
            if !Self::has_structured_syntax(*factor) {
                return Err(TensorInferenceError::invalid(
                    "chain factors must contain structured tensors",
                ));
            }
            self.validate_builtin_placeholder_channels(*factor, input, output)?;
            interfaces.push(factor_interface(self, *factor)?);
        }
        let spectators = Self::merge_explicit_interface_sequence(&interfaces)?;
        let boundary = SymbolicTensor::composite_interface(
            exposed.into_iter().chain(spectators.logical_slots()),
        );
        // A closed chain binds its endpoints locally. Its powers copy
        // the surviving boundary, without multiplying that internal pair.
        Self::merge_explicit_interface_sequence(&[boundary])
    }

    fn infer_view(&mut self, atom: AtomView<'_>) -> InferenceResult<PartialStructure> {
        #[cfg(test)]
        {
            #[cfg(debug_assertions)]
            let debug = DebugInference::is_active();
            #[cfg(not(debug_assertions))]
            let debug = false;
            if debug {
                tests::DEBUG_INFERENCE_VISITS.with(|count| count.set(count.get() + 1));
            } else {
                tests::INFERENCE_VISITS.with(|count| count.set(count.get() + 1));
            }
        }
        let interface = self.infer_view_uncached(atom)?;
        // Merging checks each regional boundary even when the caller only
        // consumes the independent verdict on the root.
        let boundary = Self::merge_explicit_interface_sequence(slice::from_ref(&interface))?;
        if self.purpose == ObservationPurpose::Validate {
            return Ok(interface);
        }
        if matches!(atom, AtomView::Fun(function)
            if composition::is_composite_head(function.get_symbol(), &SPENSO_TAG)
                || function.get_symbol() == SPENSO_TAG.dot)
        {
            self.admitted_leaf_interfaces.insert(
                atom.to_owned(),
                (
                    true,
                    boundary
                        .logical_slots()
                        .into_iter()
                        .map(composition::port_atom)
                        .collect(),
                ),
            );
        }
        self.admitted_interfaces.insert(atom, boundary);
        Ok(interface)
    }

    fn infer_view_uncached(&mut self, atom: AtomView<'_>) -> InferenceResult<PartialStructure> {
        #[cfg(feature = "reference-cases")]
        let _phase = timing::scope(Phase::InterfaceDiscovery);
        if self.leaf_inference == LeafInference::ObserveStorage {
            return self.infer_storage_view(atom);
        }
        if let Some(interface) = self.cached_interface(atom) {
            return Ok(interface.clone());
        }

        if let AtomView::Add(sum) = atom {
            let mut expected = None;
            let mut has_scalar_term = false;
            for term in sum.iter() {
                if !Self::has_structured_syntax(term) {
                    if self.purpose == ObservationPurpose::Publish {
                        self.admitted_interfaces
                            .insert(term, PartialStructure::from_logical_slots([]));
                        self.admitted_leaf_interfaces
                            .insert(term.to_owned(), (false, Vec::new()));
                    }
                    has_scalar_term = true;
                    continue;
                }
                let actual = Self::merge_explicit_interface_sequence(&[self.infer_view(term)?])?;
                let Some(current) = &expected else {
                    expected = Some(actual);
                    continue;
                };
                if !Self::additive_interfaces_match(current, &actual)
                    || (self.leaf_inference == LeafInference::ObserveEncoded
                        && !Self::encoded_open_orders_match(current, &actual))
                {
                    return Err(TensorInferenceError::invalid(
                        "tensor summands do not have compatible tensor interfaces",
                    ));
                }
            }
            let Some(expected) = expected else {
                return Err(TensorInferenceError::invalid(
                    "expression does not contain valid tagged Spenso tensor syntax",
                ));
            };
            if has_scalar_term && !expected.canonical().is_scalar() {
                return Err(TensorInferenceError::invalid(
                    "cannot add a scalar expression to a non-scalar tensor",
                ));
            }
            return Ok(expected.into_algebra_result());
        }

        if let AtomView::Mul(product) = atom {
            let mut interfaces = Vec::new();
            for factor in product.iter() {
                if self.cached_interface(factor).is_some() || Self::has_structured_syntax(factor) {
                    interfaces.push(self.infer_view(factor)?);
                } else if self.purpose == ObservationPurpose::Publish {
                    self.admitted_interfaces
                        .insert(factor, PartialStructure::from_logical_slots([]));
                    self.admitted_leaf_interfaces
                        .insert(factor.to_owned(), (false, Vec::new()));
                }
            }
            if interfaces.is_empty() {
                return Err(TensorInferenceError::invalid(
                    "expression does not contain valid tagged Spenso tensor syntax",
                ));
            }
            return Self::merge_explicit_interface_sequence(&interfaces)
                .map(TensorInterface::into_algebra_result);
        }

        if let AtomView::Pow(power) = atom {
            let (base, exponent) = power.get_base_exp();
            if !Self::has_structured_syntax(base) {
                return Err(TensorInferenceError::invalid(
                    "expression does not contain valid tagged Spenso tensor syntax",
                ));
            }
            let interface = self.infer_view(base)?;
            if Self::has_structured_syntax(exponent)
                && !self.infer_view(exponent)?.canonical().is_scalar()
            {
                return Err(TensorInferenceError::invalid(
                    "a tensor exponent must be scalar",
                ));
            }
            return Self::power_interface(&interface, exponent, || atom.to_string());
        }

        if let AtomView::Fun(function) = atom {
            let symbol = function.get_symbol();
            let arguments = function.iter().collect::<Vec<_>>();
            let unmaterialized_tensor_leaf = Self::is_unmaterialized_tensor_leaf(atom);
            if symbol == *SYM || symbol == *ANTISYM || symbol == *CYCLIC {
                let mut interfaces = Vec::new();
                for argument in arguments {
                    if !Self::has_structured_syntax(argument) {
                        continue;
                    }
                    interfaces.push(self.infer_view(argument)?);
                }
                if interfaces.is_empty() {
                    return Err(TensorInferenceError::invalid(
                        "tensor projectors require at least one tensor operand",
                    ));
                }
                return Self::merge_explicit_interface_sequence(&interfaces)
                    .map(TensorInterface::into_algebra_result);
            }

            if symbol == SPENSO_TAG.bracket {
                let mut interfaces = Vec::new();
                for argument in arguments {
                    if !Self::has_structured_syntax(argument) {
                        continue;
                    }
                    interfaces.push(self.infer_view(argument)?);
                }
                // A nonempty ordered product may consist entirely of scalar
                // factors. The empty product was rejected by syntax validation.
                return Self::merge_explicit_interface_sequence(&interfaces)
                    .map(TensorInterface::into_algebra_result);
            }

            let compact_metric = symbol == ETS.metric
                && arguments.len() == 2
                && arguments.iter().all(|argument| {
                    self.slots
                        .parse::<LibraryRep, AbstractIndex>(*argument)
                        .is_err()
                        && self
                            .slots
                            .parse_representation::<LibraryRep>(*argument)
                            .is_err()
                        && !Self::is_chain_placeholder(*argument)
                        && Self::has_structured_syntax(*argument)
                });
            if symbol == SPENSO_TAG.dot || compact_metric {
                let [left, right] = arguments.as_slice() else {
                    return Err(TensorInferenceError::invalid(
                        "inner products require exactly two operands",
                    ));
                };
                if !Self::has_structured_syntax(*left) || !Self::has_structured_syntax(*right) {
                    return Err(TensorInferenceError::invalid(
                        "inner products require rank-one tensor operands",
                    ));
                }
                let compact_representations = self
                    .slots
                    .compact_inner_product_parts::<AbstractIndex>(function)
                    .map(|(_, _, left, right)| [left, right]);
                let mut operand_interface = |argument| {
                    if let AtomView::Fun(vector) = argument
                        && self.compact_vector_port(argument)
                    {
                        // Like a gamma's consumed vector, a dot operand has no
                        // component to materialize. Observe its declared port;
                        // constructing a fresh component would replay its hook.
                        self.validate_tensor_metadata(vector)?;
                        if let Some(interface) =
                            self.direct_leaf_interface(vector, LeafInference::Materialize)
                        {
                            // Keep the existing callback-free reusable admission.
                            // Consumed callback-bearing ports remain observation-only.
                            self.cache_interface(argument, &interface);
                            return Ok(interface);
                        }
                        self.direct_leaf_interface(vector, LeafInference::Observe)
                            .ok_or_else(|| {
                                TensorInferenceError::invalid(
                                    "inner-product vector has no observable port",
                                )
                            })
                    } else {
                        self.infer_view(argument)
                    }
                };
                let left_interface = operand_interface(*left)?;
                let right_interface = operand_interface(*right)?;
                if self.purpose == ObservationPurpose::Publish {
                    self.admitted_interfaces
                        .insert(*left, left_interface.clone());
                    self.admitted_interfaces
                        .insert(*right, right_interface.clone());
                }
                let left = left_interface;
                let right = right_interface;
                if symbol == SPENSO_TAG.dot
                    && (left.logical_slots().len() != 1 || right.logical_slots().len() != 1)
                {
                    return Err(TensorInferenceError::invalid(
                        "dot requires rank-one operands; use a metric or explicit contraction for partial contractions",
                    ));
                }
                let reusable = arguments
                    .iter()
                    .all(|argument| self.reusable_interfaces.contains_key(argument.get_data()));
                let (left, left_spectators) = Self::consume_inner_product_axis(
                    &left,
                    compact_representations.map(|reps| reps[0]),
                )?;
                let (right, right_spectators) = Self::consume_inner_product_axis(
                    &right,
                    compact_representations.map(|reps| reps[1]),
                )?;
                if !left.rep().matches(&right.rep()) {
                    return Err(TensorInferenceError::invalid(
                        "inner-product operands carry incompatible representations",
                    ));
                }
                if matches!(
                    (left.aind, right.aind),
                    (PartialIndex::Explicit(left), PartialIndex::Explicit(right)) if left != right
                ) {
                    return Err(TensorInferenceError::invalid(
                        "inner-product operands carry unequal explicit indices",
                    ));
                }
                let interface =
                    Self::merge_explicit_interface_sequence(&[left_spectators, right_spectators])?;
                if reusable {
                    // A dot consumes only its paired axes. Cached operands establish that
                    // repeating their inference has no materialization effects.
                    self.cache_interface(atom, &interface);
                }
                return Ok(interface);
            }

            if symbol.has_tag(&SPENSO_TAG.broadcast) {
                let argument = arguments[0];
                if !Self::has_structured_syntax(argument) {
                    return Err(TensorInferenceError::invalid(format!(
                        "broadcast function `{symbol}` does not contain a structured tensor argument"
                    )));
                }
                return self
                    .infer_view(argument)
                    .map(TensorInterface::into_algebra_result);
            }

            if let Some(structure) = self.builtin_tensor_structure(symbol, &arguments, false)? {
                return Ok(self.builtin_interface(atom, &arguments, &structure));
            }

            if symbol == SPENSO_TAG.chain {
                return self.chain_interface_with(&arguments, Self::infer_view);
            } else if symbol == SPENSO_TAG.trace {
                let representation = self
                    .slots
                    .parse_representation::<LibraryRep>(arguments[0])
                    .map_err(|_| {
                        TensorInferenceError::invalid(
                            "trace metadata is not a Spenso representation",
                        )
                    })?;
                let dual = representation.dual();
                let mut interfaces = Vec::new();
                for factor in shadowing::trace_factor_views(&arguments[1..]) {
                    if !Self::has_structured_syntax(factor) {
                        return Err(TensorInferenceError::invalid(
                            "trace factors must contain structured tensors",
                        ));
                    }
                    self.validate_builtin_placeholder_channels(factor, representation, dual)?;
                    interfaces.push(self.infer_view(factor)?);
                }
                return Self::merge_explicit_interface_sequence(&interfaces)
                    .map(TensorInterface::into_algebra_result);
            } else if symbol.has_tag(&SPENSO_TAG.tensor) {
                self.validate_tensor_metadata(function)?;
                if unmaterialized_tensor_leaf {
                    let _ = self.direct_leaf_interface(function, self.leaf_inference);
                    return Ok(PartialStructure::from_logical_slots([]));
                }
            }
        }

        if !Self::has_structured_syntax(atom) {
            return Err(TensorInferenceError::invalid(
                "expression does not contain valid tagged Spenso tensor syntax",
            ));
        }

        if let AtomView::Fun(function) = atom
            && let Some(interface) = self.direct_leaf_interface(function, self.leaf_inference)
        {
            // OpenPortIds are local positions, canonicalized when interfaces are
            // combined. Direct parsing allocates no dummy or callback-sensitive state.
            self.cache_interface(atom, &interface);
            return Ok(interface);
        }

        self.infer_storage_view(atom)
    }

    // The existing structural parser owns stored leaf axes, including trace
    // spectator incidences that the logical interface contracts locally.
    fn infer_storage_view(&mut self, atom: AtomView<'_>) -> InferenceResult<PartialStructure> {
        let mut open_markers = HashMap::new();
        // Constructors retain their historical compact-port materialization.
        // A result check observes the already normalized expression instead:
        // replaying a leaf callback with a fresh dummy could change its rank.
        let materialized = if self.leaf_inference != LeafInference::Materialize {
            AtomOrView::View(atom)
        } else {
            AtomOrView::Atom(atom.replace_map(|value, _, output| {
                if self.slots.parse::<LibraryRep, AbstractIndex>(value).is_ok() {
                    return;
                }
                if let Ok(representation) = self.slots.parse_representation::<LibraryRep>(value) {
                    let marker = self.dummies.fresh_index();
                    open_markers.insert(marker, value.to_owned());
                    **output = representation.slot::<AbstractIndex, _>(marker).to_atom();
                }
            }))
        };
        let infer = |slots: &mut SlotMatcher| {
            OrderedStructure::<LibraryRep, AbstractIndex>::structure_from_atom(
                materialized.as_view(),
                slots,
            )
        };
        let (inferred, written) = if self.purpose == ObservationPurpose::Publish {
            self.slots.observe_slots(infer)
        } else {
            (infer(&mut self.slots), Vec::new())
        };
        let inferred = match inferred {
            Ok(inferred) => inferred,
            Err(StructureError::EmptyStructure(_)) => {
                return Ok(PartialStructure::from_logical_slots([]));
            }
            Err(error) => {
                return Err(TensorInferenceError::invalid(format!(
                    "invalid Spenso expression: {error}"
                )));
            }
        };
        // Spenso retains source order alongside canonical storage, including
        // structural bundles. No second walk of the leaf is needed.
        let logical = inferred
            .layout()
            .canonical_to_logical(&inferred.canonical().external_structure());
        if self.purpose == ObservationPurpose::Publish {
            let written = written
                .into_iter()
                .filter_map(|atom| {
                    self.slots
                        .parse::<LibraryRep, AbstractIndex>(atom.as_view())
                        .ok()
                        .map(|slot| {
                            let original = open_markers.get(&slot.aind()).cloned().unwrap_or(atom);
                            (slot, original)
                        })
                })
                .collect::<HashMap<_, _>>();
            if let Some(slots) = logical
                .iter()
                .map(|slot| written.get(slot).cloned())
                .collect::<Option<Vec<_>>>()
            {
                self.admitted_leaf_interfaces
                    .insert(atom.to_owned(), (true, slots));
            }
        }
        Ok(PartialStructure::from_logical_slots(
            logical.into_iter().map(|slot| {
                let index = if open_markers.contains_key(&slot.aind()) {
                    PartialIndex::open(0)
                } else {
                    Self::partial_index(slot.aind(), self.leaf_inference)
                };
                slot.rep().slot(index)
            }),
        ))
    }
}

impl InterfaceInference {
    fn unambiguous_anonymous_ports(interface: &PartialStructure) -> bool {
        let mut unresolved = HashSet::new();
        interface.logical_slots().into_iter().all(|slot| {
            !matches!(
                slot.aind,
                PartialIndex::Open(_) | PartialIndex::Explicit(AbstractIndex::Open { .. })
            ) || unresolved.insert(slot.rep())
        })
    }

    fn encoded_open_orders_match(left: &PartialStructure, right: &PartialStructure) -> bool {
        let open_order = |interface: &PartialStructure| {
            let mut order = HashMap::<_, Vec<_>>::new();
            for slot in interface.logical_slots() {
                if let PartialIndex::Explicit(index @ AbstractIndex::Open { .. }) = slot.aind {
                    order.entry(slot.rep()).or_default().push(index);
                }
            }
            order
        };
        open_order(left) == open_order(right)
    }

    pub fn additive_interfaces_match(left: &PartialStructure, right: &PartialStructure) -> bool {
        // Unresolved ports remain positional; only explicit indices identify ports across terms.
        let left_slots = left.logical_slots();
        let mut right_slots = right.logical_slots();
        if left_slots == right_slots {
            return true;
        }
        if left_slots.len() != right_slots.len()
            || left_slots
                .iter()
                .chain(&right_slots)
                .any(|slot| matches!(slot.aind, PartialIndex::Open(_)))
        {
            return false;
        }
        if left.canonical() == right.canonical() {
            return true;
        }
        // Symbolic dimensions need not have a total canonical ordering. Compare
        // complete explicit ports by equality, retaining multiplicity and the
        // caller's logical order instead of trusting their storage permutation.
        for slot in left_slots {
            let Some(position) = right_slots.iter().position(|candidate| *candidate == slot) else {
                return false;
            };
            right_slots.swap_remove(position);
        }
        true
    }

    pub fn has_structured_syntax(value: AtomView<'_>) -> bool {
        value.is_tensorial(StrictTensorFilter::Tagged)
            || Self::is_unmaterialized_tensor_leaf(value)
            || matches!(
                value,
                AtomView::Fun(projector)
                    if projector.get_symbol() == *shadowing::SYM
                        || projector.get_symbol() == *shadowing::ANTISYM
                        || projector.get_symbol() == *shadowing::CYCLIC
            )
    }

    /// Scan all syntax, including scalar metadata and exponents. Tensor-free
    /// input is opaque at admission; other tags and callbacks do not exclude it.
    fn contains_tensor_symbol(value: AtomView<'_>) -> bool {
        let mut found = false;
        value.visitor(&mut |node| {
            if !found {
                found = match node {
                    AtomView::Var(variable) => variable.get_symbol().has_tag(&SPENSO_TAG.tensor),
                    AtomView::Fun(function) => {
                        let symbol = function.get_symbol();
                        symbol.has_tag(&SPENSO_TAG.tensor)
                            || symbol == SPENSO_TAG.chain
                            || symbol == SPENSO_TAG.trace
                    }
                    _ => false,
                };
            }
            !found
        });
        found
    }

    fn is_unmaterialized_tensor_leaf(value: AtomView<'_>) -> bool {
        let AtomView::Fun(function) = value else {
            return false;
        };
        function.get_symbol().has_tag(&SPENSO_TAG.tensor)
            && function.iter().any(Self::is_chain_placeholder)
            && function.iter().all(|argument| {
                Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_err()
                    && Representation::<LibraryRep>::try_from(argument).is_err()
            })
    }

    fn is_chain_placeholder(value: AtomView<'_>) -> bool {
        matches!(
            value,
            AtomView::Var(variable)
                if variable.get_symbol() == SPENSO_TAG.chain_in
                    || variable.get_symbol() == SPENSO_TAG.chain_out
        )
    }

    /// Namespaced `in`/`out` placeholders are wiring labels, not free abstract indices. They may
    /// only occur as direct ports of tensor leaves inside a chain or trace factor.
    pub fn validate_placeholder_scope(
        value: AtomView<'_>,
        inside_factor: bool,
    ) -> Result<(), String> {
        Self::observe_indices(
            value,
            &mut SlotMatcher::default(),
            &mut AHashMap::new(),
            ObservationScope {
                inside_factor,
                validate: true,
                ..ObservationScope::default()
            },
            &mut |_| false,
        )
        .map(|_| ())
        .map_err(|error| error.to_string())
    }

    fn admitted_oriented_square(
        proofs: &TensorProofs,
        value: AtomView<'_>,
        base: &PartialStructure,
        exponent: AtomView<'_>,
    ) -> Option<PartialStructure> {
        if exponent != Atom::num(2).as_view()
            || base
                .logical_slots()
                .iter()
                .all(|slot| slot.rep().rep.is_self_dual())
            || !base.open_positions().is_empty()
        {
            return None;
        }
        let prior = proofs.retained_interface(value)?;
        let copies = Self::merge_explicit_interface_sequence(&[base.clone(), base.clone()]).ok()?;
        Self::additive_interfaces_match(&prior, &copies).then_some(copies)
    }

    pub(crate) fn power_interface(
        interface: &PartialStructure,
        exponent: AtomView<'_>,
        description: impl Fn() -> String,
    ) -> InferenceResult<PartialStructure> {
        if interface.canonical().is_scalar() {
            return Ok(interface.clone());
        }
        if !interface
            .logical_slots()
            .iter()
            .all(|slot| slot.rep().rep.is_self_dual())
        {
            return Err(TensorInferenceError::invalid(format!(
                "invalid power of non-self-dual tensor `{}`",
                description()
            )));
        }
        let exponent = Rational::try_from(exponent).map_err(|_| {
            TensorInferenceError::invalid(format!("invalid tensor power `{}`", description()))
        })?;
        if exponent.denominator() != 1 {
            return Err(TensorInferenceError::invalid(format!(
                "fractional tensor power `{}` has no well-defined interface",
                description()
            )));
        }
        let repetitions = exponent.numerator().abs();
        let slots = interface.logical_slots();
        for slot in &slots {
            let PartialIndex::Explicit(index) = slot.aind else {
                continue;
            };
            let compatible = slots
            .iter()
            .filter(|candidate| {
                matches!(candidate.aind, PartialIndex::Explicit(candidate_index) if candidate_index == index)
                && slot.rep().matches(&candidate.rep())
            })
            .count();
            let exceeds_einstein_multiplicity = match repetitions.to_i64() {
                Some(0) => false,
                Some(1) => compatible > 2,
                Some(2) => compatible > 1,
                _ => compatible > 0,
            };
            if exceeds_einstein_multiplicity {
                return Err(TensorInferenceError::invalid(format!(
                    "explicit index `{index}` occurs on more than two compatible ports in `{}`; \
                     to raise a contraction, bracket it first, as in `bracket(Q(i)*Q(i))^2`, \
                     `bracket(T(i,i))^-2` or `dot(Q,Q)^-2`",
                    description()
                )));
            }
        }
        if exponent.numerator() % 2 == 0 {
            return Ok(PartialStructure::from_logical_slots([]));
        }
        Ok(interface.clone().into_algebra_result())
    }

    /// Symbolica distributes an integer power over a product, so `(Q(i)*R(i))^-2`
    /// and `Q(i)^-2*R(i)^-2` are one Atom. A compatible explicit index carried by
    /// two inverse tensor powers of one product therefore has no unique pairing.
    /// Each inverse alone closes its own pair, so other factors may reuse it.
    fn reject_shared_inverse_power_indices<'a>(
        factors: impl Iterator<Item = AtomView<'a>>,
    ) -> InferenceResult<()> {
        let mut inverse_slots: Vec<PartialSlot> = Vec::new();
        for factor in factors {
            let AtomView::Pow(power) = factor else {
                continue;
            };
            let (base, exponent) = power.get_base_exp();
            let inverse = Rational::try_from(exponent)
                .is_ok_and(|exponent| exponent.numerator().is_negative());
            if !inverse || !Self::has_structured_syntax(base) {
                continue;
            }
            let base_interface = InterfaceInference::default().infer_validated(base)?;
            let interface = Self::merge_explicit_interface_sequence(&[base_interface])?;
            for slot in interface.logical_slots() {
                let PartialIndex::Explicit(index) = slot.aind else {
                    continue;
                };
                if inverse_slots.iter().any(|other| {
                    matches!(other.aind, PartialIndex::Explicit(other_index) if other_index == index)
                        && slot.rep().matches(&other.rep())
                }) {
                    return Err(TensorInferenceError::invalid(format!(
                        "explicit index `{index}` occurs in several inverse tensor powers of one \
                         product, which is also how `(Q(i)*R(i))^-2` reads after normalization; \
                         write `bracket(Q(i)*R(i))^-2` to invert a contraction, or give each \
                         inverse its own index, as in `Q(i)^-2*R(j)^-2`"
                    )));
                }
                inverse_slots.push(slot);
            }
        }
        Ok(())
    }

    /// Delay replacement storage until a child changes; the common explicit
    /// tensor leaves have no powers to lower. Every child is still visited once,
    /// preserving callback execution and error order.
    fn lower_tensor_power_arguments<'a>(
        arguments: impl Iterator<Item = AtomView<'a>> + Clone,
    ) -> InferenceResult<Option<Vec<AtomOrView<'a>>>> {
        let original = arguments.clone();
        let mut arguments = arguments.enumerate();
        while let Some((position, argument)) = arguments.next() {
            if let Some(lowered) = Self::lower_tensor_powers(argument)? {
                let mut replacements = original
                    .take(position)
                    .map(AtomOrView::View)
                    .collect::<Vec<_>>();
                replacements.push(AtomOrView::Atom(lowered));
                for (_, argument) in arguments {
                    replacements.push(match Self::lower_tensor_powers(argument)? {
                        Some(lowered) => AtomOrView::Atom(lowered),
                        None => AtomOrView::View(argument),
                    });
                }
                return Ok(Some(replacements));
            }
        }
        Ok(None)
    }

    /// Lower tensor powers, returning no replacement when the normalized atom is unchanged.
    pub fn lower_tensor_powers(value: AtomView<'_>) -> InferenceResult<Option<Atom>> {
        let lowered = match value {
            AtomView::Add(sum) => {
                let Some(terms) = Self::lower_tensor_power_arguments(sum.iter())? else {
                    return Ok(None);
                };
                terms.iter().map(AtomCore::as_atom_view).sum()
            }
            AtomView::Mul(product) => {
                Self::reject_shared_inverse_power_indices(product.iter())?;
                let Some(factors) = Self::lower_tensor_power_arguments(product.iter())? else {
                    return Ok(None);
                };
                factors.iter().map(AtomCore::as_atom_view).product()
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let lowered_base = Self::lower_tensor_powers(base)?;
                let lowered_exponent = Self::lower_tensor_powers(exponent)?;
                let base = lowered_base.as_ref().map_or(base, Atom::as_view);
                let exponent = lowered_exponent.as_ref().map_or(exponent, Atom::as_view);
                if !Self::has_structured_syntax(base) {
                    if lowered_base.is_none() && lowered_exponent.is_none() {
                        return Ok(None);
                    }
                    base.pow(exponent)
                } else {
                    // Validate multiplicity before lowering through tensor-aware multiplication.
                    let mut inference = InterfaceInference {
                        purpose: ObservationPurpose::Publish,
                        ..InterfaceInference::default()
                    };
                    let (raw_interface, occurrences) = inference.infer_observed(base)?;
                    if Self::has_structured_syntax(exponent)
                        && !inference.infer_view(exponent)?.canonical().is_scalar()
                    {
                        return Err(TensorInferenceError::invalid(
                            "a tensor exponent must be scalar",
                        ));
                    }
                    Self::power_interface(&raw_interface, exponent, || value.to_string())?;
                    let interface = Self::merge_explicit_interface_sequence(&[raw_interface])?;
                    if interface.canonical().is_scalar() {
                        base.pow(exponent)
                    } else {
                        let exponent = Rational::try_from(exponent).map_err(|_| {
                            TensorInferenceError::invalid(
                                "a non-scalar tensor power requires an integer exponent",
                            )
                        })?;
                        // A negative even power is the reciprocal of the contracted
                        // positive power, which `power_interface` has already found
                        // to be scalar. Odd negative powers keep an open interface
                        // and have no reciprocal.
                        let negative = exponent.numerator().is_negative();
                        if exponent.denominator() != 1
                            || (negative && exponent.numerator() % 2 != 0)
                        {
                            return Err(TensorInferenceError::invalid(
                                "a non-scalar tensor power requires a non-negative or even integer exponent",
                            ));
                        }
                        let repetitions =
                            usize::try_from(exponent.numerator().abs()).map_err(|_| {
                                TensorInferenceError::invalid(
                                    "tensor power exponent is too large to materialize",
                                )
                            })?;
                        if negative {
                            // The typed power owns the reciprocal, so the raw and
                            // typed spellings admit the same closed inverse.
                            let exponent = SymbolicTensor::<PartialStructure>::from_produced_parts(
                                Atom::num(exponent),
                                PartialStructure::from_logical_slots([]),
                            )?;
                            SymbolicTensor::<PartialStructure>::from_produced_parts(
                                base.to_owned(),
                                interface,
                            )?
                            .pow(&exponent)?
                            .expression
                        } else if repetitions == 0 {
                            Atom::num(1)
                        } else if interface.open_positions().is_empty() {
                            // Explicit labels already identify the contractions; only open
                            // ports need the pair selection of tensor multiplication.
                            base.pow(Atom::num(repetitions))
                        } else {
                            let factor = inference.publish_admission(
                                base.to_owned(),
                                interface,
                                occurrences,
                            );
                            let factor = factor.finish_parts()?;
                            let mut result = factor.clone();
                            for _ in 1..repetitions {
                                result = result.multiply(&factor)?;
                            }
                            result.expression
                        }
                    }
                }
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if SPENSO_TAG.is_scalar_metadata(symbol) {
                    return Ok(None);
                }
                let arguments = Self::lower_tensor_power_arguments(function.iter())?;
                if arguments.is_none()
                    && symbol.get_normalization_function().is_none()
                    && symbol.get_evaluation_info().is_none()
                {
                    return Ok(None);
                }
                // Preserve existing callback invocation when rebuilding is observable.
                let builder = FunctionBuilder::new(symbol);
                match arguments {
                    Some(arguments) => builder.add_args(arguments).finish(),
                    None => builder.add_args(function.iter()).finish(),
                }
            }
            _ => return Ok(None),
        };
        Ok((lowered.as_view() != value).then_some(lowered))
    }

    /// Select the axis consumed by a compact inner product and retain each
    /// spectator in the producer's logical order.
    pub(super) fn consume_inner_product_axis(
        interface: &PartialStructure,
        representation: Option<Representation<LibraryRep>>,
    ) -> InferenceResult<(PartialSlot, PartialStructure)> {
        let mut slots = interface.logical_slots();
        let mut compact = slots.iter().enumerate().filter_map(|(position, slot)| {
            matches!(
                slot.aind,
                PartialIndex::Open(_) | PartialIndex::Explicit(AbstractIndex::Open { .. })
            )
            .then_some(position)
        });
        let position = match (compact.next(), compact.next()) {
            (Some(position), None) => position,
            (None, None) if slots.len() == 1 => 0,
            _ => {
                return Err(TensorInferenceError::invalid(
                    "inner-product operands require one compact axis or one explicit port",
                ));
            }
        };
        let port = slots.remove(position);
        if representation.is_some_and(|rep| rep != port.rep()) {
            return Err(TensorInferenceError::invalid(
                "inner-product compact axis disagrees with its observed interface",
            ));
        }
        Ok((port, PartialStructure::from_logical_slots(slots)))
    }

    pub fn merge_explicit_interface_sequence(
        interfaces: &[PartialStructure],
    ) -> InferenceResult<PartialStructure> {
        let slots = interfaces
            .iter()
            .flat_map(PartialStructureExt::logical_slots)
            .collect::<Vec<_>>();
        let mut contracted = HashSet::new();
        for (position, slot) in slots.iter().enumerate() {
            let PartialIndex::Explicit(index) = slot.aind else {
                continue;
            };
            let compatible = slots
            .iter()
            .enumerate()
            .filter(|(candidate_position, candidate)| {
                if *candidate_position == position {
                    return true;
                }
                matches!(candidate.aind, PartialIndex::Explicit(candidate_index) if candidate_index == index)
                    && slot.rep().matches(&candidate.rep())
            })
            .map(|(candidate_position, _)| candidate_position)
            .collect::<Vec<_>>();
            if compatible.len() > 2 {
                return Err(TensorCompositionError::InvalidExplicitMultiplicity {
                    index,
                    occurrences: compatible.len(),
                }
                .into());
            }
            if compatible.len() == 2 {
                contracted.extend(compatible);
            }
        }

        // The root admission pass must retain an untouched leaf's declaration.
        // Enclosing algebra can consume that existing layout without creating
        // another temporary pair of sorting permutations.
        if interfaces.len() == 1 && contracted.is_empty() {
            return Ok(interfaces[0].clone());
        }
        let survivors = slots
            .into_iter()
            .enumerate()
            .filter(|(position, _)| !contracted.contains(position))
            .map(|(_, slot)| slot);
        Ok(SymbolicTensor::composite_interface(survivors))
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::shorthands::bracket::BracketNormalizer;
    use crate::shorthands::schoonschip::{Schoonschip, SchoonschipSettings};
    use crate::tensor::SymbolicNetParse;
    use crate::{CookSettings, test_support};
    use spenso::network::library::FunctionLibrary;
    use spenso::network::library::function_lib::Wrap;
    use spenso::network::parsing::ParseSettings;
    use spenso::structure::representation::ExtendibleReps;
    use spenso::structure::slot::DualSlotTo;
    use spenso::tensor_symbol;
    use std::cell::Cell;
    use std::mem;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Arc, Mutex};
    use symbolica::coefficient::CoefficientView;
    use symbolica::symbol;

    mod authoritative_metadata;
    mod auto_axis;
    mod builtin_words;
    mod polynomial_expansion;
    mod replacement_observer;
    mod scalar_admission;

    thread_local! {
        pub(in crate::tensor) static OCCURRENCE_FUNCTION_VISITS: Cell<usize> = const { Cell::new(0) };
        pub(super) static INTRINSIC_HEAD_CHECKS: Cell<usize> = const { Cell::new(0) };
        pub(crate) static INFERENCE_CALLS: Cell<usize> = const { Cell::new(0) };
        pub(crate) static INFERENCE_VISITS: Cell<usize> = const { Cell::new(0) };
        pub(crate) static DEBUG_INFERENCE_VISITS: Cell<usize> = const { Cell::new(0) };
        pub(crate) static DEBUG_INTRINSIC_HEAD_CHECKS: Cell<usize> = const { Cell::new(0) };
        pub(crate) static DEBUG_LEAF_PROOF_CALLS: Cell<usize> = const { Cell::new(0) };
        pub(crate) static DEBUG_INFERENCE_CALLS: Cell<usize> = const { Cell::new(0) };
        pub(crate) static DEBUG_SCOPE_VALIDATIONS: Cell<usize> = const { Cell::new(0) };
        pub(crate) static DEBUG_OCCURRENCE_FUNCTION_VISITS: Cell<usize> = const { Cell::new(0) };
        pub(super) static LEAF_PROOF_CALLS: Cell<usize> = const { Cell::new(0) };
        pub(crate) static SCOPE_VALIDATIONS: Cell<usize> = const { Cell::new(0) };
    }

    #[test]
    fn producer_results_keep_callback_ports_and_history_without_release_inference() {
        test_support::test_initialize();
        let representation = LibraryRep::from(ColorFundamental {}).new_rep(3);
        let slot = representation.slot::<AbstractIndex, _>(AbstractIndex::Normal(99711));
        let head = tensor_symbol!("producer_transport::callback", norm = |_node, _out| {});
        let source =
            SymbolicTensor::infer(FunctionBuilder::new(head).add_arg(slot.to_atom()).finish())
                .unwrap();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        SCOPE_VALIDATIONS.with(|count| count.set(0));
        OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
        DEBUG_INFERENCE_CALLS.with(|count| count.set(0));
        for finish in [
            SymbolicTensor::<PartialStructure>::with_produced_expression,
            SymbolicTensor::<PartialStructure>::with_algebra_result,
            SymbolicTensor::<PartialStructure>::with_rewritten_expression,
        ] {
            let result = finish(&source, Atom::num(2) * source.expression()).unwrap();
            assert_eq!(result.structure, source.structure);
            assert!(
                result
                    .explicit_index_occurrences()
                    .unwrap()
                    .is_bounded_by(source.explicit_index_occurrences().unwrap())
            );
            assert!(
                source
                    .explicit_index_occurrences()
                    .unwrap()
                    .is_bounded_by(result.explicit_index_occurrences().unwrap())
            );
            let zero = finish(&result, Atom::Zero).unwrap();
            assert_eq!(zero.rank(), 1);
            assert_eq!(zero.structure, source.structure);
            assert!(Arc::ptr_eq(
                zero.proofs.occurrences.get().unwrap(),
                result.proofs.occurrences.get().unwrap(),
            ));
        }
        let identity = source
            .with_identity_result(Atom::num(3) * source.expression(), None)
            .unwrap();
        assert_eq!(identity.structure, source.structure);
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
        assert!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()) > 0);
        #[cfg(debug_assertions)]
        assert!(DEBUG_INFERENCE_CALLS.with(|count| count.get()) > 0);
        assert!(source.with_checked_expression(Atom::one()).is_err());
    }

    #[test]
    fn producer_slot_maps_transport_erased_oriented_history() {
        test_support::test_initialize();
        let representation = LibraryRep::from(ColorFundamental {}).new_rep(3);
        let free = representation.slot::<AbstractIndex, _>(AbstractIndex::Normal(99721));
        let hidden = representation.slot::<AbstractIndex, _>(AbstractIndex::Dummy(99723));
        let head = tensor_symbol!("producer_transport::oriented");
        let leaf = |slot: Slot<LibraryRep, AbstractIndex>| {
            FunctionBuilder::new(head).add_arg(slot.to_atom()).finish()
        };
        let source = SymbolicTensor::infer(
            leaf(free) * leaf(hidden) * leaf(hidden.rep().dual().slot(hidden.aind())),
        )
        .unwrap();
        let zero = source.with_produced_expression(Atom::Zero).unwrap();
        let dual = PartialStructure::from_logical_slots([free
            .rep()
            .dual()
            .slot(PartialIndex::Explicit(free.aind()))]);
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
        let result = zero
            .with_produced_parts(Atom::Zero, dual.clone(), |slot| {
                slot.rep().dual().slot(slot.aind())
            })
            .unwrap();
        assert_eq!(result.structure, dual);
        assert_eq!(
            result.explicit_index_occurrences().unwrap().count(hidden),
            2
        );
        assert_eq!(
            result
                .explicit_index_occurrences()
                .unwrap()
                .count(free.rep().dual().slot(free.aind())),
            1
        );
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn materializer_accounts_new_dummies_without_release_interface_inference() {
        test_support::test_initialize();
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(3);
        let free = representation.slot::<AbstractIndex, _>(AbstractIndex::Normal(99731));
        let dummy = representation.slot::<AbstractIndex, _>(AbstractIndex::Dummy(99733));
        let left = tensor_symbol!("producer_transport::materialized_left");
        let right = tensor_symbol!("producer_transport::materialized_right");
        let source =
            SymbolicTensor::infer(FunctionBuilder::new(left).add_arg(free.to_atom()).finish())
                .unwrap();
        let expression = source.expression()
            * FunctionBuilder::new(left).add_arg(dummy.to_atom()).finish()
            * FunctionBuilder::new(right)
                .add_arg(dummy.to_atom())
                .finish();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let result = source
            .with_materialized_expression("test materializer", expression)
            .unwrap();
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.explicit_index_occurrences().unwrap().count(dummy), 2);
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn callback_output_requires_an_admitted_leaf_or_valid_tensor_syntax() {
        test_support::test_initialize();
        let port = ExtendibleReps::EUCLIDEAN
            .new_rep(3)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(99749))
            .to_atom();
        let head = tensor_symbol!("producer_transport::callback_grammar", norm = |_, _| {});
        let source = SymbolicTensor::infer(symbolica::function!(head, &port, &port)).unwrap();
        assert!(source.is_scalar());
        let invalid = source
            .with_produced_expression(Atom::var(head))
            .unwrap_err();
        assert!(
            invalid
                .to_string()
                .contains("must be called with at least one structural port")
        );
    }

    #[test]
    fn callback_dummy_renaming_accounts_history_without_interface_inference() {
        test_support::test_initialize();
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(3);
        let free = representation.slot::<AbstractIndex, _>(AbstractIndex::Normal(99751));
        let before = representation.slot::<AbstractIndex, _>(AbstractIndex::Dummy(99753));
        let after = representation.slot::<AbstractIndex, _>(AbstractIndex::Dummy(99755));
        let old = before.to_atom();
        let new = after.to_atom();
        let head = tensor_symbol!(
            "producer_transport::rename_dummy",
            norm = move |node, out| {
                if let AtomView::Fun(function) = node {
                    let arguments = function.iter().collect::<Vec<_>>();
                    if arguments
                        .first()
                        .is_some_and(|value| *value == Atom::num(1).as_view())
                        && arguments
                            .get(2)
                            .is_some_and(|value| *value == old.as_view())
                    {
                        **out = FunctionBuilder::new(function.get_symbol())
                            .add_arg(arguments[0])
                            .add_arg(arguments[1])
                            .add_arg(&new)
                            .add_arg(&new)
                            .finish();
                    }
                }
            }
        );
        let source = SymbolicTensor::infer(
            FunctionBuilder::new(head)
                .add_arg(0)
                .add_arg(free.to_atom())
                .add_arg(before.to_atom())
                .add_arg(before.to_atom())
                .finish(),
        )
        .unwrap();
        let rewritten = FunctionBuilder::new(head)
            .add_arg(1)
            .add_arg(free.to_atom())
            .add_arg(before.to_atom())
            .add_arg(before.to_atom())
            .finish();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let result = source.with_produced_expression(rewritten).unwrap();
        assert_eq!(result.structure, source.structure);
        assert_eq!(
            result.explicit_index_occurrences().unwrap().count(before),
            2
        );
        assert_eq!(result.explicit_index_occurrences().unwrap().count(after), 2);
        let leaves = result.proofs.leaf_interfaces.get().unwrap().lock().unwrap();
        assert_eq!(
            leaves[&result.expression].1,
            [free.to_atom(), after.to_atom(), after.to_atom()]
        );
        assert!(!leaves.contains_key(&source.expression));
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "tensor producer introduced unreported index occurrences")]
    fn intrinsic_producer_must_report_allocated_dummies() {
        test_support::test_initialize();
        let source = SymbolicTensor::infer(Atom::one()).unwrap();
        let dummy = ExtendibleReps::EUCLIDEAN
            .new_rep(3)
            .slot::<AbstractIndex, _>(AbstractIndex::Dummy(99761));
        let left = FunctionBuilder::new(tensor_symbol!("producer_transport::unreported_left"))
            .add_arg(dummy.to_atom())
            .finish();
        let right = FunctionBuilder::new(tensor_symbol!("producer_transport::unreported_right"))
            .add_arg(dummy.to_atom())
            .finish();
        source.with_produced_expression(left * right).unwrap();
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "tensor producer contract violated")]
    fn producer_interface_mismatch_is_a_debug_assertion() {
        test_support::test_initialize();
        let slot = ExtendibleReps::EUCLIDEAN
            .new_rep(3)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(99741));
        let source = SymbolicTensor::infer(
            FunctionBuilder::new(tensor_symbol!("producer_transport::invalid"))
                .add_arg(slot.to_atom())
                .finish(),
        )
        .unwrap();
        source.with_produced_expression(Atom::one()).unwrap();
    }

    #[test]
    fn callback_identity_results_discard_terminal_certificates_for_changed_payloads() {
        test_support::test_initialize();
        let scalar_callback =
            symbol!("producer_transport::closed_work"; Scalar; norm = |_node, _out| {});
        let source =
            SymbolicTensor::infer(FunctionBuilder::new(scalar_callback).add_arg(0).finish())
                .unwrap();
        source.reduction_observations();
        let index = ExtendibleReps::MINKOWSKI
            .new_rep(4)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(99771));
        let trace = shadowing::trace(
            Bispinor {}.new_rep(4).to_symbolic([]),
            [
                FunctionBuilder::new(AGS.gamma)
                    .add_args([
                        Atom::var(SPENSO_TAG.chain_in),
                        Atom::var(SPENSO_TAG.chain_out),
                        index.to_atom(),
                    ])
                    .finish(),
                FunctionBuilder::new(AGS.gamma)
                    .add_args([
                        Atom::var(SPENSO_TAG.chain_in),
                        Atom::var(SPENSO_TAG.chain_out),
                        index.to_atom(),
                    ])
                    .finish(),
            ],
        );
        let result = source
            .with_identity_result(
                trace,
                Some(DomainObservations::terminal_polynomial(
                    &source.structure,
                    false,
                )),
            )
            .unwrap();
        assert_eq!(result.structure, source.structure);
        assert!(!result.reduction_observations().is_terminal_polynomial());
        assert_ne!(result.reduction_observations().identity_candidates(), 0);
    }

    #[test]
    fn certified_identity_does_not_certify_new_callback_heads() {
        crate::test_support::test_initialize();
        let source = SymbolicTensor::infer(Atom::one()).unwrap();
        let callback = symbolica::symbol!("trusted_identity_new_callback"; Scalar;
            norm = |_value, _out| {});
        let result = source
            .with_identity_result(FunctionBuilder::new(callback).add_arg(7).finish(), None)
            .unwrap();
        assert!(result.proofs.is_admitted());
        assert!(!result.normalization_is_intrinsic());
        let tensor = FunctionBuilder::new(spenso::tensor_symbol!("trusted_identity_added_port"))
            .add_arg(spenso::mink!(4, 98491))
            .finish();
        assert!(result.with_checked_expression(tensor).is_err());
    }

    #[test]
    fn shared_expansion_preserves_variable_polynomial_dispatch_and_identity() {
        let x = Atom::var(symbolica::symbol!("shared_expansion_x"));
        let y = Atom::var(symbolica::symbol!("shared_expansion_y"));
        let f = FunctionBuilder::new(symbolica::symbol!("shared_expansion_f"))
            .add_arg(&y)
            .finish();
        let source = (&x + &f).pow(3) + (&x + &y) * (&f - &y);
        let typed = SymbolicTensor::<PartialStructure>::from_produced_parts(
            source.clone(),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        for var in [None, Some(x.as_view()), Some(f.as_view())] {
            for via_poly in [false, true] {
                let expected = match (var, via_poly) {
                    (Some(var), true) => source.expand_via_poly::<i32, Atom>(Some(var.to_owned())),
                    (Some(var), false) => source.expand_in(var),
                    (None, true) => source.expand_via_poly::<i32, Atom>(None),
                    (None, false) => source.expand(),
                };
                let result = typed.expanded(var, via_poly).unwrap();
                assert_eq!(result.expression, expected);
                assert_eq!(result.structure, typed.structure);
            }
        }
        for via_poly in [false, true] {
            assert_eq!(
                typed
                    .expanded(Some(Atom::num(2).as_view()), via_poly)
                    .unwrap_err()
                    .to_string(),
                "Expansion must be done wrt an indeterminate"
            );
        }
        let mut unchanged = SymbolicTensor::<PartialStructure>::from_produced_parts(
            source.expand(),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        unchanged.is_metric = true;
        unchanged.is_composite = false;
        assert_eq!(unchanged.expanded(None, false).unwrap(), unchanged);
        assert_eq!(typed.expression, source);
    }

    #[test]
    fn polynomial_expansion_accepts_exponents_larger_than_i16() {
        let x = Atom::var(symbolica::symbol!("large_expansion_x"));
        let y = Atom::var(symbolica::symbol!("large_expansion_y"));
        let source = x.pow(32768) * (&x + &y);
        let expected = x.pow(32769) + x.pow(32768) * &y;
        let typed = SymbolicTensor::<PartialStructure>::infer(source).unwrap();
        for var in [None, Some(x.as_view())] {
            let result = typed.expanded(var, true).unwrap();
            assert_eq!(result.expression, expected);
            assert_eq!(result.structure, typed.structure);
        }
    }

    #[test]
    fn shared_expansion_preserves_tensor_ports_zero_and_callback_fallback() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let slots = [821, 823].map(|i| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(i))
                .to_atom()
        });
        let leaf = FunctionBuilder::new(spenso::tensor_symbol!("shared_expand_tensor"))
            .add_args(&slots)
            .finish();
        let x = Atom::var(symbolica::symbol!("shared_expand_port_x"));
        let y = Atom::var(symbolica::symbol!("shared_expand_port_y"));
        for source in [
            (&x + &y) * &leaf + &x * &leaf,
            (&x + &y) * &leaf - &x * &leaf - &y * &leaf,
        ] {
            let mut typed = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
            typed.structure = PartialStructure::from_logical_slots(
                typed.structure.logical_slots().into_iter().rev(),
            );
            let result = typed.expanded(None, false).unwrap();
            assert_eq!(result.expression, source.expand());
            assert_eq!(result.structure, typed.structure);
            assert_eq!(typed.expression, source);
        }
        let unresolved = FunctionBuilder::new(spenso::tensor_symbol!("shared_expand_unresolved"))
            .add_arg(rep.to_symbolic([]))
            .finish();
        let source = (&x + &y) * unresolved;
        let typed = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
        assert!(!typed.structure.open_positions().is_empty());
        let expected = typed.with_algebra_result(source.expand()).unwrap();
        assert_eq!(typed.expanded(None, false).unwrap(), expected);

        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let hook = symbolica::symbol!(
            "shared_expand_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let scalar = FunctionBuilder::new(symbolica::symbol!("shared_expand_opaque"; Scalar))
            .add_arg(FunctionBuilder::new(hook).add_arg(&leaf).finish())
            .finish();
        let source = (&x + &y) * scalar + Atom::num(1);
        let typed = SymbolicTensor::<PartialStructure>::from_produced_parts(
            source.clone(),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        calls.store(0, Ordering::Relaxed);
        let expected = typed.with_algebra_result(source.expand()).unwrap();
        let old_calls = calls.load(Ordering::Relaxed);
        calls.store(0, Ordering::Relaxed);
        assert_eq!(typed.expanded(None, false).unwrap(), expected);
        assert_eq!(calls.load(Ordering::Relaxed), old_calls);
    }

    #[test]
    fn expansion_keeps_closed_dummy_pairs_between_their_factors() {
        test_support::test_initialize();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let [mu, nu] = [99763, 99767].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let heads = [
            spenso::vector_symbol!("closed_expand_p"),
            spenso::vector_symbol!("closed_expand_q"),
        ];
        let vector =
            |head: Symbol, index: &Atom| FunctionBuilder::new(head).add_arg(index).finish();
        let p = |index| vector(heads[0], index);
        let q = |index| vector(heads[1], index);
        let x = Atom::var(symbol!("closed_expand_x"));
        let y = Atom::var(symbol!("closed_expand_y"));
        let dot = |index| p(index) * q(index);
        let canonical =
            |value: &SymbolicTensor<PartialStructure>| value.canonize().unwrap().expression;
        let closed = |index| spenso::bracket!(dot(index));
        // Each source closes `p·q` in a sum, while the same name also occurs
        // in a neighbouring factor or in the copies of an expanded power.
        let cases = [
            (
                (&x + dot(&mu)).pow(2),
                &x * &x + Atom::num(2) * &x * dot(&mu) + closed(&mu).pow(2),
            ),
            ((&x + dot(&mu)) * p(&mu), &x * p(&mu) + p(&mu) * dot(&nu)),
            (
                &y + (&x + dot(&mu)).pow(2) * dot(&mu),
                &y + &x * &x * dot(&mu)
                    + Atom::num(2) * &x * dot(&mu) * dot(&nu)
                    + dot(&mu) * closed(&nu).pow(2),
            ),
        ];
        for (source, expected) in cases {
            let typed = SymbolicTensor::infer(source).unwrap();
            let expected = SymbolicTensor::infer(expected).unwrap();
            for (var, via_poly) in [
                (None, false),
                (None, true),
                (Some(x.as_view()), false),
                (Some(x.as_view()), true),
            ] {
                let result = typed.expanded(var, via_poly).unwrap();
                assert_eq!(result.structure, typed.structure);
                // The published payload must itself be admissible with the
                // same contractions; renamed dummies compare canonically.
                let admitted = SymbolicTensor::infer(result.expression.clone()).unwrap();
                assert_eq!(admitted.structure, typed.structure);
                let complete = result.expanded(None, false).unwrap();
                let complete = SymbolicTensor::infer(complete.expression.clone()).unwrap();
                assert_eq!(canonical(&complete), canonical(&expected));
            }
        }

        // Expansion can also cancel the sum which closed a pair, leaving it in
        // the outer product scope of the result.
        let typed = SymbolicTensor::infer((&x + dot(&mu)) * &y - &x * &y).unwrap();
        let exposed = typed.expanded(None, false).unwrap();
        assert_eq!(exposed.expression, &y * dot(&mu));
        let probe = SymbolicTensor::infer(p(&mu)).unwrap();
        assert!(typed.multiply(&probe).is_ok());
        assert!(matches!(
            exposed.multiply(&probe),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));
    }

    /// A closed `p·q` written with one dummy name.
    fn closed_dot(index: usize, heads: [Symbol; 2]) -> Atom {
        let slot = ExtendibleReps::MINKOWSKI
            .new_rep(Dimension::Concrete(4))
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom();
        heads
            .map(|head| FunctionBuilder::new(head).add_arg(&slot).finish())
            .iter()
            .fold(Atom::one(), |product, vector| product * vector)
    }

    fn summands(value: &SymbolicTensor<PartialStructure>) -> usize {
        match value.expression.as_view() {
            AtomView::Add(sum) => sum.get_nargs(),
            _ => 1,
        }
    }

    fn assert_unchanged(
        result: &SymbolicTensor<PartialStructure>,
        source: &SymbolicTensor<PartialStructure>,
    ) {
        assert_eq!(result.expression, source.expression);
        assert_eq!(result.structure, source.structure);
        assert_eq!(
            result
                .explicit_index_occurrences()
                .unwrap()
                .entries()
                .collect::<Vec<_>>(),
            source
                .explicit_index_occurrences()
                .unwrap()
                .entries()
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn function_selector_marking_keeps_opaque_scalar_callbacks_untouched() {
        test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let hook = symbol!("scope_function_opaque_hook"; Scalar; norm=move |_,_| {
            observed.fetch_add(1, Ordering::Relaxed);
        });
        let q = FunctionBuilder::new(spenso::vector_symbol!("scope_opaque_q"))
            .add_arg(spenso::mink!(4, 99923))
            .finish();
        let raw = FunctionBuilder::new(hook).add_arg(&q).finish();
        let source = SymbolicTensor::infer(raw.clone()).unwrap();
        for via_poly in [false, true] {
            calls.store(0, Ordering::Relaxed);
            let result = source.expanded(Some(q.as_view()), via_poly).unwrap();
            assert_eq!(result.expression, source.expression);
            assert_eq!(
                calls.load(Ordering::Relaxed),
                0,
                "selected occurrences inside a backend-opaque scalar function are not visited"
            );
        }
    }

    #[test]
    fn function_selector_carrier_does_not_reserve_user_symbols() {
        test_support::test_initialize();
        let p = FunctionBuilder::new(spenso::vector_symbol!("scope_carrier_p"))
            .add_arg(spenso::mink!(4, 99925))
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("scope_carrier_q"))
            .add_arg(spenso::mink!(4, 99925))
            .finish();
        let r = FunctionBuilder::new(spenso::vector_symbol!("scope_carrier_r"))
            .add_arg(spenso::mink!(4, 99925))
            .finish();
        let user = symbol!("expansion_selector");
        let x = Atom::var(symbol!("scope_carrier_x"));
        let source = SymbolicTensor::infer(Atom::var(user) * (&x + p * &q) * r).unwrap();
        for via_poly in [false, true] {
            let result = source.expanded(Some(q.as_view()), via_poly).unwrap();
            assert!(result.expression.contains_symbol(user));
            assert_eq!(summands(&result), 2);
        }
    }

    #[test]
    fn separated_function_selector_does_not_visit_opaque_arguments() {
        test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let hook = symbol!("scope_separated_opaque_hook"; Scalar; norm=move |_,_| {
            observed.fetch_add(1, Ordering::Relaxed);
        });
        let mu = spenso::mink!(4, 99931);
        let p = FunctionBuilder::new(spenso::vector_symbol!("scope_sep_p"))
            .add_arg(&mu)
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("scope_sep_q"))
            .add_arg(&mu)
            .finish();
        let r = FunctionBuilder::new(spenso::vector_symbol!("scope_sep_r"))
            .add_arg(&mu)
            .finish();
        let x = Atom::var(symbol!("scope_separated_x"));
        let opaque = FunctionBuilder::new(hook).add_arg(&q).finish();
        let source = SymbolicTensor::infer((&x + p * &q) * &r + opaque * r).unwrap();
        for via_poly in [false, true] {
            calls.store(0, Ordering::Relaxed);
            let result = source.expanded(Some(q.as_view()), via_poly).unwrap();
            assert_eq!(result.structure, source.structure);
            assert_eq!(
                calls.load(Ordering::Relaxed),
                0,
                "a separate scope rename does not open an opaque function"
            );
        }
    }

    #[test]
    fn separated_function_selector_does_not_repeat_renaming_callbacks() {
        test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let q_head = spenso::vector_symbol!(
            "scope_repeated_q",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let mu = spenso::mink!(4, 99933);
        let p = FunctionBuilder::new(spenso::vector_symbol!("scope_repeat_p"))
            .add_arg(&mu)
            .finish();
        let q = FunctionBuilder::new(q_head).add_arg(&mu).finish();
        let r = FunctionBuilder::new(spenso::vector_symbol!("scope_repeat_r"))
            .add_arg(&mu)
            .finish();
        let x = Atom::var(symbol!("scope_repeated_x"));
        let source = SymbolicTensor::infer((&x + p * &q) * r).unwrap();
        for via_poly in [false, true] {
            calls.store(0, Ordering::Relaxed);
            let scalar_selected = source.expanded(Some(x.as_view()), via_poly).unwrap();
            let once = calls.load(Ordering::Relaxed);
            calls.store(0, Ordering::Relaxed);
            let function_selected = source.expanded(Some(q.as_view()), via_poly).unwrap();
            let actual = calls.load(Ordering::Relaxed);
            assert_eq!(function_selected.structure, scalar_selected.structure);
            assert_eq!(
                actual, once,
                "marking must not first perform and discard a rename"
            );
        }
    }

    #[test]
    fn function_selector_has_no_fixed_candidate_name_budget() {
        test_support::test_initialize();
        let occupied = (0..1024)
            .map(|index| {
                Atom::var(
                    SymbolBuilder::new(NamespacedSymbol::parse(&format!(
                        "idenso::expansion_selector_{index}"
                    )))
                    .build()
                    .unwrap(),
                )
            })
            .collect::<Vec<_>>();
        let coefficient = Atom::mul_many(&occupied);
        let mu = spenso::mink!(4, 99935);
        let nu = spenso::mink!(4, 99937);
        let p_head = spenso::vector_symbol!("scope_budget_p");
        let q_head = spenso::vector_symbol!("scope_budget_q");
        let vector =
            |head: Symbol, index: &Atom| FunctionBuilder::new(head).add_arg(index).finish();
        let p = vector(p_head, &mu);
        let q = vector(q_head, &mu);
        let r = FunctionBuilder::new(spenso::vector_symbol!("scope_budget_r"))
            .add_arg(&mu)
            .finish();
        let x = Atom::var(symbol!("scope_budget_x"));
        let source = SymbolicTensor::infer(&coefficient * (&x + p * &q) * &r).unwrap();
        assert!(
            source.expanded(Some(q.as_view()), false).is_ok(),
            "private carrier selection must not reserve a finite name set from users"
        );
        // The carrier skips every occupied name, and removing it keeps each
        // user symbol and the separated value.
        let carrier = SymbolicTensor::<PartialStructure>::expansion_carrier(
            source.expression.as_view(),
            q.as_view(),
        )
        .unwrap();
        assert!(
            occupied
                .iter()
                .all(|symbol| !symbol.contains_symbol(carrier))
        );
        let expected =
            &coefficient * &x * &r + &coefficient * vector(p_head, &nu) * vector(q_head, &nu) * &r;
        let canonical = |value: &Atom| {
            SymbolicTensor::infer(value.clone())
                .unwrap()
                .canonize()
                .unwrap()
                .expression
        };
        for via_poly in [false, true] {
            let expanded = source.expanded(Some(q.as_view()), via_poly).unwrap();
            assert!(!expanded.expression.contains_symbol(carrier));
            for symbol in &occupied {
                assert!(
                    expanded.expression.contains(symbol.as_view()),
                    "poly={via_poly}"
                );
            }
            assert_eq!(summands(&expanded), 2);
            assert_eq!(canonical(&expanded.expression), canonical(&expected));
        }
    }

    #[test]
    fn function_selector_inside_opaque_leaves_still_selects_enclosing_sums() {
        test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let hook = symbol!("scope_parent_opaque_hook"; Scalar; norm=move |_,_| {
            observed.fetch_add(1, Ordering::Relaxed);
        });
        let mu = spenso::mink!(4, 99941);
        let p = FunctionBuilder::new(spenso::vector_symbol!("scope_parent_p"))
            .add_arg(&mu)
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("scope_parent_q"))
            .add_arg(&mu)
            .finish();
        let r = FunctionBuilder::new(spenso::vector_symbol!("scope_parent_r"))
            .add_arg(&mu)
            .finish();
        let x = Atom::var(symbol!("scope_parent_x"));
        let y = Atom::var(symbol!("scope_parent_y"));
        let opaque = FunctionBuilder::new(hook).add_arg(&q).finish();
        let source = SymbolicTensor::infer((&x + p * &q) * r * (&y + &opaque)).unwrap();
        for via_poly in [false, true] {
            calls.store(0, Ordering::Relaxed);
            let result = source.expanded(Some(q.as_view()), via_poly).unwrap();
            assert_eq!(
                summands(&result),
                if via_poly { 2 } else { 4 },
                "poly={via_poly} output={}",
                result.expression
            );
            assert_eq!(result.structure, source.structure);
            assert_eq!(calls.load(Ordering::Relaxed), 0);
            assert!(result.expression.contains(opaque.as_view()));
            SymbolicTensor::infer(result.expression).unwrap();
        }
    }

    #[test]
    fn selected_closed_scalar_powers_remove_carriers_inside_brackets() {
        test_support::test_initialize();
        let mu = spenso::mink!(4, 99943);
        let p = FunctionBuilder::new(spenso::vector_symbol!("scope_power_p"))
            .add_arg(&mu)
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("scope_power_q"))
            .add_arg(&mu)
            .finish();
        let x = Atom::var(symbol!("scope_power_x"));
        let source = SymbolicTensor::infer((&x + p * &q).pow(12)).unwrap();
        let carrier = SymbolicTensor::<PartialStructure>::expansion_carrier(
            source.expression.as_view(),
            q.as_view(),
        )
        .unwrap();
        for via_poly in [false, true] {
            let result = source.expanded(Some(q.as_view()), via_poly).unwrap();
            assert_eq!(
                summands(&result),
                13,
                "poly={via_poly} output={}",
                result.expression
            );
            assert!(!result.expression.contains_symbol(carrier));
            let expected = source.expanded(Some(x.as_view()), via_poly).unwrap();
            assert_eq!(
                result.canonize().unwrap().expression,
                expected.canonize().unwrap().expression
            );
            SymbolicTensor::infer(result.expression).unwrap();
        }
    }

    #[test]
    fn selected_scalar_coefficients_of_closed_products_remain_visible() {
        test_support::test_initialize();
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(1));
        let [mu, nu] = [99951, 99953].map(|i| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(i))
                .to_atom()
        });
        let heads = [
            spenso::vector_symbol!("selected_coefficient_p"),
            spenso::vector_symbol!("selected_coefficient_q"),
            spenso::vector_symbol!("selected_coefficient_r"),
            spenso::vector_symbol!("selected_coefficient_s"),
        ];
        let vector = |head, index: &Atom| FunctionBuilder::new(head).add_arg(index).finish();
        let first = vector(heads[0], &mu) * vector(heads[1], &mu);
        let second = vector(heads[2], &nu) * vector(heads[3], &nu);
        let x = Atom::var(symbol!("selected_coefficient_x"));
        let y = Atom::var(symbol!("selected_coefficient_y"));
        for (raw, expected) in [
            ((&x * &first + &y).pow(4), 50625),
            ((&y + (&x + &first) * &second).pow(2), 80089),
        ] {
            let source = SymbolicTensor::infer(raw).unwrap();
            for via_poly in [false, true] {
                let result = source.expanded(Some(x.as_view()), via_poly).unwrap();
                assert!(summands(&result) > 1);
                assert_eq!(result.structure, source.structure);
                SymbolicTensor::infer(result.expression.clone()).unwrap();
                let numeric = result.expression.replace_map(|node, _, out| {
                    if node == x.as_view() {
                        **out = Atom::num(2);
                    } else if node == y.as_view() {
                        **out = Atom::num(3);
                    } else if let AtomView::Fun(fun) = node
                        && let Some(position) =
                            heads.iter().position(|head| *head == fun.get_symbol())
                    {
                        **out = Atom::num([2, 3, 5, 7][position]);
                    }
                });
                assert_eq!(
                    BracketNormalizer::normalize(numeric.as_view()),
                    Atom::num(expected),
                    "poly={via_poly}"
                );
                // The selected scalar stays outside every introduced bracket.
                result.expression.visitor(&mut |node| {
                    if let AtomView::Fun(fun) = node
                        && fun.get_symbol() == SPENSO_TAG.bracket
                    {
                        assert!(
                            !node.contains(x.as_view()),
                            "poly={via_poly} output={}",
                            result.expression
                        );
                        return false;
                    }
                    true
                });
            }
        }
    }

    #[test]
    fn selected_closed_coefficient_sums_keep_their_pairs_when_raised() {
        test_support::test_initialize();
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let mu = rep
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(99961))
            .to_atom();
        let vector = |head: Symbol| FunctionBuilder::new(head).add_arg(&mu).finish();
        let p = vector(spenso::vector_symbol!("raised_coefficient_p"));
        let q = vector(spenso::vector_symbol!("raised_coefficient_q"));
        let [x, y, z] =
            ["x", "y", "z"].map(|name| Atom::var(symbol!(&format!("raised_coefficient_{name}"))));
        // The closed sum multiplies `x` inside a raised term: its own pair
        // `p·q` must stay one contraction in each copy, never `p^2 q^2`.
        let raw = (Atom::one() + &x * (&x + &y * &p * &q + &z)).pow(Atom::num(2));
        let source = SymbolicTensor::infer(raw.clone()).unwrap();
        for via_poly in [false, true] {
            let result = source.expanded(Some(x.as_view()), via_poly).unwrap();
            SymbolicTensor::infer(result.expression.clone()).unwrap();
            result.expression.visitor(&mut |node| {
                if let AtomView::Pow(power) = node {
                    assert!(
                        !matches!(power.get_base_exp().0, AtomView::Fun(fun)
                            if fun.get_symbol() == p.get_symbol().unwrap()
                                || fun.get_symbol() == q.get_symbol().unwrap()),
                        "poly={via_poly} output={}",
                        result.expression
                    );
                }
                true
            });
            // One generated name renames the closed pair; with it fixed, the
            // result is the expansion of `(1 + x^2 + x·y·(p·q) + x·z)^2`.
            let mut generated = HashSet::new();
            result.expression.visitor(&mut |node| {
                if let AtomView::Fun(fun) = node
                    && (fun.get_symbol() == p.get_symbol().unwrap()
                        || fun.get_symbol() == q.get_symbol().unwrap())
                {
                    generated.insert(fun.iter().next().unwrap().to_owned());
                    return false;
                }
                true
            });
            assert_eq!(generated.len(), 1, "poly={via_poly}");
            let generated = generated.into_iter().next().unwrap();
            let fixed = result.expression.replace_map(|node, _, out| {
                if node == generated.as_view() {
                    **out = mu.clone();
                }
            });
            // Each closed term is one bracketed scalar `b = y·(p·q)`.
            let b = spenso::bracket!(&y * &p * &q);
            let two = Atom::num(2);
            let expected = Atom::one()
                + &two * &x * &z
                + &two * &x * &b
                + &two * x.pow(Atom::num(2))
                + &two * x.pow(Atom::num(2)) * &z * &b
                + x.pow(Atom::num(2)) * b.pow(Atom::num(2))
                + x.pow(Atom::num(2)) * z.pow(Atom::num(2))
                + &two * x.pow(Atom::num(3)) * &z
                + &two * x.pow(Atom::num(3)) * &b
                + x.pow(Atom::num(4));
            // The polynomial route keeps coefficients of `x` factored. Only a
            // powered bracket scopes its pair; drop the others to compare.
            let complete = |value: Atom| {
                SymbolicTensor::infer(value)
                    .unwrap()
                    .expanded(None, false)
                    .unwrap()
                    .expression
                    .replace_map(|node, context, out| {
                        if let AtomView::Fun(fun) = node
                            && fun.get_symbol() == SPENSO_TAG.bracket
                            && context.parent_type != Some(symbolica::atom::AtomType::Pow)
                        {
                            **out = fun.iter().next().unwrap().to_owned();
                        }
                    })
            };
            assert_eq!(complete(fixed), complete(expected), "poly={via_poly}");
        }
    }

    #[test]
    fn expansion_keeps_arguments_of_opaque_untagged_functions() {
        test_support::test_initialize();
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let [i, j] = [99_971, 99_972].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        // Untagged heads with index-shaped arguments are opaque scalars at
        // admission: their arguments are data, not indices, so expansion must
        // neither rename them nor treat the base as a closed scalar.
        let [u0, u1, u2, u4] =
            ["u0", "u1", "u2", "u4"].map(|name| symbol!(&format!("opaque_untagged_{name}")));
        let call = |head: Symbol, index: &Atom| FunctionBuilder::new(head).add_arg(index).finish();
        let x = Atom::var(symbol!("opaque_untagged_x"));
        let sources = [
            (call(u4, &i) + call(u0, &j) * call(u1, &j) * call(u2, &i)).pow(Atom::num(2)),
            (&x + call(u0, &i) * call(u1, &i)) * call(u2, &i),
            (&x + call(u0, &i) * call(u1, &i)).pow(Atom::num(3)) * call(u4, &i),
        ];
        for raw in sources {
            let source = SymbolicTensor::infer(raw.clone()).unwrap();
            assert!(source.is_scalar());
            for (var, via_poly) in [
                (None, false),
                (None, true),
                (Some(x.clone()), false),
                (Some(x.clone()), true),
                (Some(call(u4, &i)), false),
                (Some(call(u4, &i)), true),
            ] {
                let result = source
                    .expanded(var.as_ref().map(Atom::as_view), via_poly)
                    .unwrap();
                // Exactly the plain scalar expansion: every opaque call keeps
                // its literal argument.
                let expected = match &var {
                    None if via_poly => raw.expand_via_poly::<i32, Atom>(None),
                    None => raw.expand(),
                    Some(var) if via_poly => raw.expand_via_poly::<i32, Atom>(Some(var.clone())),
                    Some(var) => raw.expand_in(var.as_view()),
                };
                assert_eq!(
                    result.expression, expected,
                    "{raw} var={var:?} poly={via_poly}"
                );
                assert!(result.is_scalar());
            }
        }
    }

    #[test]
    fn wrap_indices_keeps_arguments_of_opaque_untagged_functions() {
        test_support::test_initialize();
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let [i, j] = [99_981, 99_982].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let scope = symbol!("opaque_wrap_scope");
        let call =
            |head: Symbol, argument: &Atom| FunctionBuilder::new(head).add_arg(argument).finish();
        let [u0, u1, u2, w] =
            ["u0", "u1", "u2", "w"].map(|name| symbol!(&format!("opaque_wrap_{name}")));
        let [p, q, r] = ["p", "q", "r"]
            .map(|name| SPENSO_TAG.rank_one_tensor_symbol(&format!("opaque_wrap_{name}")));
        let x = Atom::var(symbol!("opaque_wrap_x"));
        let scoped = |slot: &Atom| AbstractIndex::wrap_expression(slot.as_view(), scope, |_| true);
        let opaque = call(w, &call(p, &i));
        // Untagged heads are opaque scalars, also beside tagged tensors: their
        // arguments are data and keep their literal names. Tagged pairs are
        // still scoped.
        let cases = [
            (call(u0, &i) * call(u1, &i), call(u0, &i) * call(u1, &i)),
            (
                (&x + call(u0, &i) * call(u1, &i)) * call(u2, &i),
                (&x + call(u0, &i) * call(u1, &i)) * call(u2, &i),
            ),
            (
                &opaque * call(q, &j) * call(r, &j),
                &opaque * call(q, &scoped(&j)) * call(r, &scoped(&j)),
            ),
        ];
        for (source, expected) in cases {
            let typed = SymbolicTensor::infer(source.clone()).unwrap();
            for dummies_only in [true, false] {
                let wrapped = typed.wrap_indices(scope, dummies_only).unwrap();
                assert_eq!(
                    wrapped.expression, expected,
                    "{source} dummies_only={dummies_only}"
                );
                assert_eq!(wrapped.structure, typed.structure);
                SymbolicTensor::infer(wrapped.expression.clone()).unwrap();
            }
        }
    }

    /// Exact rational with checked `i128` parts, for component checks.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct Exact(i128, i128);

    impl Exact {
        fn new(numerator: i128, denominator: i128) -> Self {
            assert_ne!(denominator, 0, "division by zero");
            let (mut a, mut b) = (numerator.unsigned_abs(), denominator.unsigned_abs());
            while b != 0 {
                (a, b) = (b, a % b);
            }
            let divisor = a.max(1) as i128 * denominator.signum();
            Self(numerator / divisor, denominator / divisor)
        }

        fn add(self, other: Self) -> Self {
            Self::new(self.0 * other.1 + other.0 * self.1, self.1 * other.1)
        }

        fn mul(self, other: Self) -> Self {
            Self::new(self.0 * other.0, self.1 * other.1)
        }

        fn pow(self, exponent: i64) -> Self {
            let base = if exponent < 0 {
                Self::new(self.1, self.0)
            } else {
                self
            };
            (0..exponent.unsigned_abs()).fold(Self(1, 1), |power, _| power.mul(base))
        }
    }

    /// A component value that is scalar or open on one two-dimensional index.
    #[derive(Clone, Copy, Debug)]
    enum Part {
        Scalar(Exact),
        Open([Exact; 2]),
    }

    impl Part {
        /// The product, contracting two open values over their common index.
        fn mul(self, other: Self) -> Self {
            match (self, other) {
                (Self::Scalar(a), Self::Scalar(b)) => Self::Scalar(a.mul(b)),
                (Self::Scalar(a), Self::Open(b)) | (Self::Open(b), Self::Scalar(a)) => {
                    Self::Open(b.map(|entry| a.mul(entry)))
                }
                (Self::Open(a), Self::Open(b)) => Self::Scalar(a[0].mul(b[0]).add(a[1].mul(b[1]))),
            }
        }
    }

    /// Components of rank-one and rank-two Euclidean tensors in two dimensions.
    /// A vector raised to an even power `n` is `(b·b)^(n/2)`; an odd positive
    /// power is `b (b·b)^((n-1)/2)`; a product sums each name two factors
    /// export, and a leaf that writes a name twice is its own trace.
    struct VectorComponents {
        heads: Vec<Symbol>,
        matrices: Vec<Symbol>,
        scalars: Vec<Symbol>,
    }

    impl VectorComponents {
        fn component(&self, head: Symbol, index: usize) -> Exact {
            let position = self.heads.iter().position(|h| *h == head).unwrap() as i128;
            Exact::new(
                [3, -2, 5, 1, 4][position as usize] + 7 * index as i128 * (position + 1),
                1,
            )
        }

        fn matrix(&self, head: Symbol, row: usize, column: usize) -> Exact {
            let position = self.matrices.iter().position(|h| *h == head).unwrap() as i128;
            Exact::new([2 + position, 3, 5, 7 + 2 * position][2 * row + column], 1)
        }

        fn scalar(&self, symbol: Symbol) -> Exact {
            let position = self.scalars.iter().position(|s| *s == symbol).unwrap() as i128;
            Exact::new(2 + position, 1)
        }

        fn exported(value: AtomView<'_>) -> Vec<Atom> {
            match value {
                AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.dot => vec![],
                AtomView::Fun(function) if function.get_symbol() != SPENSO_TAG.bracket => {
                    let written = function.iter().collect::<Vec<_>>();
                    written
                        .iter()
                        .filter(|name| written.iter().filter(|other| other == name).count() == 1)
                        .map(|name| name.to_owned())
                        .collect()
                }
                AtomView::Fun(function) => Self::exported(function.iter().next().unwrap()),
                AtomView::Add(sum) => Self::exported(sum.iter().next().unwrap()),
                AtomView::Mul(product) => {
                    let written = product.iter().flat_map(Self::exported).collect::<Vec<_>>();
                    written
                        .iter()
                        .filter(|name| written.iter().filter(|other| other == name).count() == 1)
                        .cloned()
                        .collect()
                }
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    let odd = i64::try_from(exponent).unwrap() % 2 != 0;
                    if odd { Self::exported(base) } else { vec![] }
                }
                _ => vec![],
            }
        }

        fn value(&self, value: AtomView<'_>, bound: &HashMap<Atom, usize>) -> Exact {
            let summed = |names: &[Atom], term: &dyn Fn(&HashMap<Atom, usize>) -> Exact| {
                (0..1_usize << names.len()).fold(Exact(0, 1), |total, code| {
                    let mut inner = bound.clone();
                    for (bit, name) in names.iter().enumerate() {
                        inner.insert(name.clone(), (code >> bit) & 1);
                    }
                    total.add(term(&inner))
                })
            };
            match value {
                AtomView::Num(number) => match number.get_coeff_view() {
                    CoefficientView::Natural(numerator, denominator, 0, 1) => {
                        Exact::new(i128::from(numerator), i128::from(denominator))
                    }
                    _ => panic!("coefficient {value}"),
                },
                AtomView::Var(variable) => self.scalar(variable.get_symbol()),
                AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.bracket => {
                    self.value(function.iter().next().unwrap(), bound)
                }
                // Admission writes a contracted pair of vectors as `dot(a, b)`.
                AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.dot => {
                    let heads = function
                        .iter()
                        .map(|vector| vector.get_symbol().unwrap())
                        .collect::<Vec<_>>();
                    (0..2).fold(Exact(0, 1), |total, c| {
                        total.add(self.component(heads[0], c).mul(self.component(heads[1], c)))
                    })
                }
                AtomView::Fun(function) if self.matrices.contains(&function.get_symbol()) => {
                    let [row, column] = [0, 1].map(|k| function.iter().nth(k).unwrap().to_owned());
                    if row == column {
                        (0..2).fold(Exact(0, 1), |total, c| {
                            total.add(self.matrix(function.get_symbol(), c, c))
                        })
                    } else {
                        self.matrix(function.get_symbol(), bound[&row], bound[&column])
                    }
                }
                AtomView::Fun(function) => {
                    let argument = function.iter().next().unwrap().to_owned();
                    self.component(function.get_symbol(), bound[&argument])
                }
                AtomView::Add(sum) => sum.iter().fold(Exact(0, 1), |total, term| {
                    total.add(self.value(term, bound))
                }),
                AtomView::Mul(product) => {
                    let factors = product.iter().collect::<Vec<_>>();
                    let written = factors
                        .iter()
                        .flat_map(|f| Self::exported(*f))
                        .collect::<Vec<_>>();
                    let mut pairs = written
                        .iter()
                        .filter(|name| written.iter().filter(|other| other == name).count() == 2)
                        .cloned()
                        .collect::<Vec<_>>();
                    pairs.sort();
                    pairs.dedup();
                    summed(&pairs, &|inner| {
                        factors.iter().fold(Exact(1, 1), |total, factor| {
                            total.mul(self.value(*factor, inner))
                        })
                    })
                }
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    let exponent = i64::try_from(exponent).unwrap();
                    let open = Self::exported(base);
                    if open.is_empty() {
                        return self.value(base, bound).pow(exponent);
                    }
                    let square = summed(&open, &|inner| self.value(base, inner).pow(2));
                    if exponent % 2 == 0 {
                        square.pow(exponent / 2)
                    } else {
                        assert!(exponent > 0, "odd negative tensor power {value}");
                        self.value(base, bound).mul(square.pow((exponent - 1) / 2))
                    }
                }
            }
        }
    }

    fn negative_power_setup() -> (VectorComponents, [Atom; 2], Atom, [Atom; 2]) {
        test_support::test_initialize();
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let heads = (0..5)
            .map(|k| SPENSO_TAG.rank_one_tensor_symbol(&format!("negative_power_v{k}")))
            .collect::<Vec<_>>();
        let matrices = (0..2)
            .map(|k| SPENSO_TAG.tensor_symbol(&format!("negative_power_m{k}")))
            .collect::<Vec<_>>();
        let scalars = ["x", "y"].map(|name| symbol!(&format!("negative_power_{name}")));
        let [i, j] = [99_991, 99_992].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let compact = rep.to_symbolic([]);
        (
            VectorComponents {
                heads,
                matrices,
                scalars: scalars.to_vec(),
            },
            [i, j],
            compact,
            scalars.map(Atom::var),
        )
    }

    #[test]
    fn inverse_boundaries_of_traces_rank_two_bases_and_brackets() {
        let (components, [i, j], _, [x, _]) = negative_power_setup();
        let v = |k: usize, port: &Atom| {
            FunctionBuilder::new(components.heads[k])
                .add_arg(port)
                .finish()
        };
        let m = |k: usize, row: &Atom, column: &Atom| {
            FunctionBuilder::new(components.matrices[k])
                .add_arg(row)
                .add_arg(column)
                .finish()
        };
        let bracket = |value: Atom| {
            FunctionBuilder::new(SPENSO_TAG.bracket)
                .add_arg(&value)
                .finish()
        };
        let dot = |left: Atom, right: Atom| {
            FunctionBuilder::new(SPENSO_TAG.dot)
                .add_arg(&left)
                .add_arg(&right)
                .finish()
        };
        let [q, r, s, u] = [0, 1, 2, 3].map(|k| v(k, &i));
        let trace = m(0, &i, &i);
        let n = |value: i64| Atom::num(value);
        // Source and the value which main's Spenso executor gives for the same
        // components. A trace closes its own pair; named axes of a rank-two base
        // pair by label; `bracket` is the user's spelling of a closed scalar.
        let accepted = [
            (&x * trace.pow(n(-1)), Part::Scalar(Exact(2, 9))),
            (m(0, &i, &j).pow(n(-2)), Part::Scalar(Exact(1, 87))),
            (
                (m(0, &i, &j) + m(1, &j, &i)).pow(n(-2)),
                Part::Scalar(Exact(1, 409)),
            ),
            (
                m(0, &i, &j).pow(n(-2)) * &r,
                Part::Open([Exact(-2, 87), Exact(4, 29)]),
            ),
            (
                (&s + &u * bracket(&q * &r).pow(n(-1))).pow(n(-2)),
                Part::Scalar(Exact(6498, 4_642_045)),
            ),
            (
                (&s + &u * dot(q.clone(), r.clone()).pow(n(-1))).pow(n(-2)),
                Part::Scalar(Exact(6498, 4_642_045)),
            ),
            (bracket(&q * &r).pow(n(-2)), Part::Scalar(Exact(1, 12_996))),
            (
                bracket(q.pow(n(2))).pow(n(-1)) * &r,
                Part::Open([Exact(-2, 109), Exact(12, 109)]),
            ),
            (bracket(&q * &q).pow(n(-2)), Part::Scalar(Exact(1, 11_881))),
            (
                bracket(trace.clone()).pow(n(-1)) * &x,
                Part::Scalar(Exact(2, 9)),
            ),
            (
                bracket(trace.clone()).pow(n(-1)) * &r,
                Part::Open([Exact(-2, 9), Exact(4, 3)]),
            ),
            (bracket(q.clone()).pow(n(-2)), Part::Scalar(Exact(1, 109))),
            (
                bracket(&q * v(1, &j)).pow(n(-2)),
                Part::Scalar(Exact(1, 16_132)),
            ),
            (
                bracket(bracket(&q * &r) * v(2, &j) * v(2, &j)).pow(n(-2)),
                Part::Scalar(Exact(1, 6_386_247_396)),
            ),
        ];
        let value =
            |expression: &Atom, main: Part| match main {
                Part::Scalar(_) => {
                    Part::Scalar(components.value(expression.as_view(), &HashMap::new()))
                }
                Part::Open(_) => Part::Open([0, 1].map(|c| {
                    components.value(expression.as_view(), &HashMap::from([(i.clone(), c)]))
                })),
            };
        for (source, main) in accepted {
            let typed = SymbolicTensor::infer(source.clone())
                .unwrap_or_else(|error| panic!("{source}: {error}"));
            let readmitted = SymbolicTensor::infer(typed.expression.clone()).unwrap();
            assert_eq!(readmitted.structure, typed.structure, "{source}");
            for result in [&typed, &readmitted] {
                assert_eq!(
                    format!("{:?}", value(&result.expression, main)),
                    format!("{main:?}"),
                    "{source} -> {}",
                    result.expression
                );
            }
        }
        // The typed power of a rank-two base keeps its labels as well.
        let minus_two = SymbolicTensor::infer(n(-2)).unwrap();
        for (base, expected) in [
            (m(0, &i, &j), Exact(1, 87)),
            (m(0, &i, &j) + m(1, &j, &i), Exact(1, 409)),
        ] {
            let raised = SymbolicTensor::infer(base)
                .unwrap()
                .pow(&minus_two)
                .unwrap();
            assert_eq!(
                components.value(raised.expression.as_view(), &HashMap::new()),
                expected
            );
        }
        // A closed scalar inside a rank-one base keeps its own pair when the
        // open port becomes automatic: (S/(Q·R))^-2 = (Q·R)^2/(S·S).
        let closed = dot(q.clone(), r.clone());
        let raised = SymbolicTensor::infer(closed.pow(n(-1)) * &s)
            .unwrap()
            .pow(&minus_two)
            .unwrap();
        assert!(
            raised.expression.contains(closed.as_view()),
            "{}",
            raised.expression
        );
        assert_eq!(
            components.value(raised.expression.as_view(), &HashMap::new()),
            Exact(12_996, 701)
        );
        // Refusals name the bracket spelling of what is meant. A trace beside a
        // factor reusing its label stays refused, as at stage 10; its bracket
        // spelling above is admitted.
        let refused = [
            (bracket(q.clone()).pow(n(-1)), "even integer exponent"),
            (trace.pow(n(-2)), "bracket(T(i,i))^-2"),
            (q.pow(n(4)), "bracket(Q(i)*Q(i))^2"),
            (q.pow(n(-2)) * r.pow(n(-2)), "bracket(Q(i)*R(i))^-2"),
            (trace.pow(n(-1)) * &r, "more than two compatible"),
        ];
        for (source, message) in refused {
            let error = SymbolicTensor::infer(source.clone())
                .unwrap_err()
                .to_string();
            assert!(error.contains(message), "{source}: {error}");
        }
    }

    #[test]
    fn inverse_squares_close_their_pair_beside_factors_with_the_same_label() {
        let (components, [i, j], _, _) = negative_power_setup();
        let v = |k: usize, port: &Atom| {
            FunctionBuilder::new(components.heads[k])
                .add_arg(port)
                .finish()
        };
        let [q, r, s] = [0, 1, 2].map(|k| v(k, &i));
        let minus_two = SymbolicTensor::infer(Atom::num(-2)).unwrap();
        let scope = symbol!("negative_power_scope");
        // The base of the inverse square, the factor beside it, and the value
        // which main's Spenso executor gives for the same components.
        let accepted = [
            (
                q.clone(),
                r.clone(),
                Part::Open([Exact(-2, 109), Exact(12, 109)]),
            ),
            (q.clone(), &r * &s, Part::Scalar(Exact(302, 109))),
            (
                q.clone(),
                r.pow(Atom::num(2)),
                Part::Scalar(Exact(148, 109)),
            ),
            (
                q.clone(),
                v(1, &j).pow(Atom::num(-2)),
                Part::Scalar(Exact(1, 16_132)),
            ),
            (
                &q + &r,
                s.clone(),
                Part::Open([Exact(1, 97), Exact(26, 485)]),
            ),
        ];
        for (base, beside, main) in accepted {
            let source = base.pow(Atom::num(-2)) * &beside;
            let value = |expression: &Atom| match main {
                Part::Scalar(_) => {
                    Part::Scalar(components.value(expression.as_view(), &HashMap::new()))
                }
                Part::Open(_) => Part::Open([0, 1].map(|c| {
                    components.value(expression.as_view(), &HashMap::from([(i.clone(), c)]))
                })),
            };
            let typed = SymbolicTensor::infer(source.clone())
                .unwrap_or_else(|error| panic!("{source}: {error}"));
            // The typed composition writes the same closed inverse; only the
            // numbering of its fresh dummies may differ, so compare values.
            let composed = SymbolicTensor::infer(base)
                .unwrap()
                .pow(&minus_two)
                .unwrap()
                .multiply(&SymbolicTensor::infer(beside).unwrap())
                .unwrap();
            assert_eq!(composed.structure, typed.structure, "{source}");
            let mut results = vec![
                composed
                    .canonize()
                    .unwrap_or_else(|error| panic!("{source} -> {}: {error}", composed.expression)),
                typed.canonize().unwrap(),
                SymbolicTensor::infer(typed.expression.clone()).unwrap(),
                typed.wrap_indices(scope, true).unwrap(),
                typed.expanded(None, false).unwrap(),
                typed.expanded(None, true).unwrap(),
                typed,
            ];
            if matches!(main, Part::Scalar(_)) {
                results.push(results[6].wrap_indices(scope, false).unwrap());
            }
            // Each form keeps the inverse closed, its interface and its value.
            for result in results {
                let readmitted = SymbolicTensor::infer(result.expression.clone())
                    .unwrap_or_else(|error| panic!("{source} -> {}: {error}", result.expression));
                assert_eq!(
                    readmitted.structure, composed.structure,
                    "{source} -> {}",
                    result.expression
                );
                assert_eq!(
                    format!("{:?}", value(&result.expression)),
                    format!("{main:?}"),
                    "{source} -> {}",
                    result.expression
                );
            }
        }
        // `(Q(i)*R(i))^-2` is distributed to the same Atom as `Q(i)^-2*R(i)^-2`,
        // which main read as 1/((Q·Q)(R·R)), not 1/(Q·R)^2. `Q(i)^-2*Q(i)` is
        // the odd power `Q(i)^-1`. A positive power exposes its label.
        let refused = [
            (
                (&q * &r).pow(Atom::num(-2)),
                "several inverse tensor powers",
            ),
            (
                q.pow(Atom::num(-2)) * r.pow(Atom::num(-2)),
                "several inverse tensor powers",
            ),
            (
                (&q + &r).pow(Atom::num(-2)) * s.pow(Atom::num(-2)),
                "several inverse tensor powers",
            ),
            (q.pow(Atom::num(-2)) * &q, "integer exponent"),
            (q.pow(Atom::num(2)) * &r, "more than two compatible"),
            (
                q.pow(Atom::num(2)) * r.pow(Atom::num(2)),
                "more than two compatible",
            ),
        ];
        for (source, message) in refused {
            let error = SymbolicTensor::infer(source.clone())
                .unwrap_err()
                .to_string();
            assert!(error.contains(message), "{source}: {error}");
        }
        // Equal inverses with distinct labels merge into one power of the same
        // scalar, and `Q(i)^-2*Q(i)^2 = 1` wherever `Q·Q` is nonzero.
        let merged =
            SymbolicTensor::infer(q.pow(Atom::num(-2)) * v(0, &j).pow(Atom::num(-2))).unwrap();
        assert_eq!(
            components.value(merged.expression.as_view(), &HashMap::new()),
            Exact(1, 11_881)
        );
        assert_eq!(
            SymbolicTensor::infer(q.pow(Atom::num(-2)) * q.pow(Atom::num(2)))
                .unwrap()
                .expression,
            Atom::num(1)
        );
    }

    #[test]
    fn negative_even_powers_of_contracted_vectors_are_reciprocal_squares() {
        let (components, [i, j], compact, [x, y]) = negative_power_setup();
        let v = |k: usize, port: &Atom| {
            FunctionBuilder::new(components.heads[k])
                .add_arg(port)
                .finish()
        };
        let dot = |k: usize, l: usize| {
            (0..2).fold(Exact(0, 1), |total, c| {
                total.add(
                    components
                        .component(components.heads[k], c)
                        .mul(components.component(components.heads[l], c)),
                )
            })
        };
        let sum_square = (0..2).fold(Exact(0, 1), |total, c| {
            let entry = components.component(components.heads[0], c).add(
                components
                    .scalar(components.scalars[0])
                    .mul(components.component(components.heads[1], c)),
            );
            total.add(entry.mul(entry))
        });
        let [ex, ey] = [0, 1].map(|k| components.scalar(components.scalars[k]));
        // Source, then its value written with explicit components.
        let cases = [
            (v(0, &i).pow(Atom::num(-2)), dot(0, 0).pow(-1)),
            (
                &x * v(0, &i).pow(Atom::num(-2)) + &y,
                ex.mul(dot(0, 0).pow(-1)).add(ey),
            ),
            (
                v(0, &i).pow(Atom::num(-2)) * v(1, &j) * v(2, &j),
                dot(0, 0).pow(-1).mul(dot(1, 2)),
            ),
            (
                (v(0, &i) + &x * v(1, &i)).pow(Atom::num(-2)),
                sum_square.pow(-1),
            ),
            (v(0, &compact).pow(Atom::num(-2)), dot(0, 0).pow(-1)),
            (v(0, &compact).pow(Atom::num(-4)), dot(0, 0).pow(-2)),
            (v(0, &compact).pow(Atom::num(-6)), dot(0, 0).pow(-3)),
            // The positive square is unchanged.
            (v(0, &i).pow(Atom::num(2)), dot(0, 0)),
        ];
        for (source, expected) in cases {
            let typed = SymbolicTensor::infer(source.clone())
                .unwrap_or_else(|error| panic!("{source}: {error}"));
            assert!(typed.is_scalar(), "{source}");
            assert_eq!(
                components.value(typed.expression.as_view(), &HashMap::new()),
                expected,
                "{source} -> {}",
                typed.expression
            );
            // The admitted payload is stable under re-admission.
            assert_eq!(
                SymbolicTensor::infer(typed.expression.clone())
                    .unwrap()
                    .expression,
                typed.expression
            );
        }
        // The inverse must close its own pair even when the surrounding
        // product uses the exact same label, once or as another pair.
        let base = v(0, &i);
        let inverse = base.pow(Atom::num(-2));
        let typed_inverse = SymbolicTensor::infer(base)
            .unwrap()
            .pow(&SymbolicTensor::infer(Atom::num(-2)).unwrap())
            .unwrap();
        assert_eq!(
            SymbolicTensor::infer(inverse.clone()).unwrap().expression,
            typed_inverse.expression
        );
        let spectator = v(1, &i);
        let open = SymbolicTensor::infer(&inverse * &spectator).unwrap();
        assert_eq!(
            open.structure.logical_slots(),
            SymbolicTensor::infer(spectator)
                .unwrap()
                .structure
                .logical_slots()
        );
        let scalar = SymbolicTensor::infer(inverse * v(1, &i) * v(2, &i)).unwrap();
        assert!(scalar.is_scalar());
        assert_eq!(
            components.value(scalar.expression.as_view(), &HashMap::new()),
            dot(0, 0).pow(-1).mul(dot(1, 2))
        );
    }

    #[test]
    fn ill_defined_negative_tensor_powers_stay_refused() {
        let (components, [i, _], compact, _) = negative_power_setup();
        let v = |k: usize, port: &Atom| {
            FunctionBuilder::new(components.heads[k])
                .add_arg(port)
                .finish()
        };
        // An odd negative power keeps an open index and has no reciprocal; an
        // explicit index used by more than one copy pair is ambiguous, as for
        // the positive powers; a fractional power has no interface; a colour
        // index is not self-dual, so its square is no contraction.
        let colour = ColorFundamental {}
            .new_rep(3)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(99_993))
            .to_atom();
        let refused = [
            (v(0, &i).pow(Atom::num(-1)), "integer exponent"),
            (
                v(0, &i).pow(Atom::num(-3)),
                "more than two compatible ports",
            ),
            (
                v(0, &i).pow(Atom::num(-4)),
                "more than two compatible ports",
            ),
            (v(0, &i).pow(Atom::num(4)), "more than two compatible ports"),
            (v(0, &compact).pow(Atom::num(-1)), "integer exponent"),
            (v(0, &compact).pow(Atom::num(-3)), "integer exponent"),
            (
                v(0, &i).pow(Atom::num(-2) / Atom::num(3)),
                "fractional tensor power",
            ),
            (v(0, &colour).pow(Atom::num(-2)), "non-self-dual"),
            (v(0, &colour).pow(Atom::num(2)), "non-self-dual"),
        ];
        for (source, message) in refused {
            let error = SymbolicTensor::infer(source.clone())
                .unwrap_err()
                .to_string();
            assert!(error.contains(message), "{source}: {error}");
        }
    }

    #[test]
    fn tensor_powers_from_minus_six_to_six_match_component_forms() {
        let (components, [i, j], compact, [x, _]) = negative_power_setup();
        let v = |k: usize, port: &Atom| {
            FunctionBuilder::new(components.heads[k])
                .add_arg(port)
                .finish()
        };
        let vector =
            |k: usize| Part::Open([0, 1].map(|c| components.component(components.heads[k], c)));
        let mut state = 0x5eed_u64;
        let mut next = |bound: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % bound
        };
        // Checked values by exponent sign, and inverse squares beside a factor
        // that reuses their label.
        let mut checked = [0, 0];
        let mut reused = 0;
        for _ in 0..800 {
            let k = next(5) as usize;
            let l = (k + 1 + next(4) as usize) % 5;
            // Two further heads for factors outside the power.
            let others = (0..5).filter(|h| *h != k && *h != l).collect::<Vec<_>>();
            let m = others[next(3) as usize];
            let n = others.iter().copied().find(|h| *h != m).unwrap();
            let exponent = next(13) as i64 - 6;
            // The base and one component of it: a vector with an explicit
            // index, a sum of two, a product with a scalar, or compact ones.
            let (base, explicit, entry): (Atom, bool, Box<dyn Fn(usize) -> Exact>) = match next(5) {
                0 => (
                    v(k, &i),
                    true,
                    Box::new(|c| components.component(components.heads[k], c)),
                ),
                1 => (
                    v(k, &i) + &x * v(l, &i),
                    true,
                    Box::new(|c| {
                        components.component(components.heads[k], c).add(
                            components
                                .scalar(components.scalars[0])
                                .mul(components.component(components.heads[l], c)),
                        )
                    }),
                ),
                2 => (
                    &x * v(k, &i),
                    true,
                    Box::new(|c| {
                        components
                            .scalar(components.scalars[0])
                            .mul(components.component(components.heads[k], c))
                    }),
                ),
                3 => (
                    &x * v(k, &compact),
                    false,
                    Box::new(|c| {
                        components
                            .scalar(components.scalars[0])
                            .mul(components.component(components.heads[k], c))
                    }),
                ),
                _ => (
                    v(k, &compact),
                    false,
                    Box::new(|c| components.component(components.heads[k], c)),
                ),
            };
            // A factor outside the power and its value. Kinds 3 to 5 reuse
            // the label `i` of an explicit base.
            let spectator_kind = next(6);
            let (spectator, spectator_part) = match spectator_kind {
                0 => (Atom::num(1), Part::Scalar(Exact(1, 1))),
                1 => (
                    x.clone(),
                    Part::Scalar(components.value(x.as_view(), &HashMap::new())),
                ),
                2 => (v(k, &j) * v(l, &j), vector(k).mul(vector(l))),
                3 => (v(m, &i), vector(m)),
                4 => (v(m, &i) * v(n, &i), vector(m).mul(vector(n))),
                _ => (
                    v(m, &i).pow(Atom::num(-2)),
                    match vector(m).mul(vector(m)) {
                        Part::Scalar(square) => Part::Scalar(square.pow(-1)),
                        Part::Open(_) => unreachable!(),
                    },
                ),
            };
            let source = base.pow(Atom::num(exponent)) * &spectator;
            let base_accepted = if explicit {
                matches!(exponent, 0 | 1 | 2 | -2)
            } else {
                exponent >= 0 || exponent % 2 == 0
            };
            // An inverse square closes its own pair, so a factor outside may
            // reuse its label. A positive power exposes its label; a second
            // inverse sharing it is ambiguous with an inverse power of a product.
            let accepted = base_accepted
                && !(explicit
                    && matches!(
                        (exponent, spectator_kind),
                        (-2, 5) | (1, 4) | (2, 3) | (2, 4)
                    ));
            let typed = SymbolicTensor::infer(source.clone());
            assert_eq!(typed.is_ok(), accepted, "{source}: {typed:?}");
            // The typed power of the admitted base follows the same rule.
            let raised = SymbolicTensor::infer(base.clone())
                .unwrap()
                .pow(&SymbolicTensor::infer(Atom::num(exponent)).unwrap());
            assert_eq!(
                raised.is_ok(),
                base_accepted,
                "{base}^{exponent}: {raised:?}"
            );
            let (Ok(typed), Ok(raised)) = (typed, raised) else {
                continue;
            };
            if !explicit && exponent % 2 != 0 {
                // Odd compact powers keep an automatic port; their admission
                // is checked above, their components by the explicit cases.
                continue;
            }
            // The component form: the base contracted with itself.
            let square = (0..2).fold(Exact(0, 1), |total, c| total.add(entry(c).mul(entry(c))));
            let power = if exponent % 2 == 0 {
                Part::Scalar(square.pow(exponent / 2))
            } else {
                Part::Open([entry(0), entry(1)])
            };
            let expected = power.mul(spectator_part);
            for c in 0..2 {
                let bound = HashMap::from([(i.clone(), c)]);
                let (actual, wanted) = match expected {
                    Part::Scalar(value) => (
                        components.value(typed.expression.as_view(), &HashMap::new()),
                        value,
                    ),
                    Part::Open(entries) => (
                        components.value(typed.expression.as_view(), &bound),
                        entries[c],
                    ),
                };
                assert_eq!(actual, wanted, "{source} -> {}", typed.expression);
            }
            if let Part::Scalar(value) = power {
                assert_eq!(
                    components.value(raised.expression.as_view(), &HashMap::new()),
                    value,
                    "{base}^{exponent} -> {}",
                    raised.expression
                );
                assert!(raised.is_scalar());
            }
            checked[usize::from(exponent < 0)] += 1;
            reused += usize::from(explicit && exponent == -2 && matches!(spectator_kind, 3 | 4));
        }
        assert!(
            checked[0] > 50 && checked[1] > 50 && reused > 10,
            "{checked:?} {reused}"
        );
    }

    #[test]
    fn tensor_sum_powers_keep_the_internal_pairs_of_each_term() {
        test_support::test_initialize();
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(1));
        let mu = rep
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(99955))
            .to_atom();
        let heads = [
            spenso::vector_symbol!("scope_open_p"),
            spenso::vector_symbol!("scope_open_q"),
            spenso::vector_symbol!("scope_open_r"),
            spenso::vector_symbol!("scope_open_v"),
        ];
        let vector = |head, index: &Atom| FunctionBuilder::new(head).add_arg(index).finish();
        for internal in [rep, LibraryRep::from(ColorFundamental {}).new_rep(1)] {
            let nu = internal.slot::<AbstractIndex, _>(AbstractIndex::Normal(99957));
            let raw = (vector(heads[3], &mu)
                + vector(heads[0], &nu.to_atom())
                    * vector(heads[1], &nu.dual().to_atom())
                    * vector(heads[2], &mu))
            .pow(2);
            let source = SymbolicTensor::infer(raw).unwrap();
            assert!(source.is_scalar());
            let selected = vector(heads[1], &nu.dual().to_atom());
            for var in [None, Some(selected.as_view())] {
                for via_poly in [false, true] {
                    let result = source.expanded(var, via_poly).unwrap();
                    SymbolicTensor::infer(result.expression.clone()).unwrap();
                    assert_eq!(
                        result.expanded(var, via_poly).unwrap().expression,
                        result.expression
                    );
                    result.expression.visitor(&mut |node| {
                        if let AtomView::Fun(fun) = node
                            && [heads[2], heads[3]].contains(&fun.get_symbol())
                        {
                            assert_eq!(fun.iter().next().unwrap(), mu.as_view());
                        }
                        true
                    });
                    let canonical = result.canonize().unwrap();
                    SymbolicTensor::infer(canonical.expression.clone()).unwrap();
                    assert_eq!(
                        canonical.canonize().unwrap().expression,
                        canonical.expression
                    );
                    let numeric = result.expression.replace_map(|node, _, out| {
                        if let AtomView::Fun(fun) = node
                            && let Some(position) =
                                heads.iter().position(|head| *head == fun.get_symbol())
                        {
                            **out = Atom::num([2, 3, 5, 7][position]);
                        }
                    });
                    assert_eq!(
                        BracketNormalizer::normalize(numeric.as_view()),
                        Atom::num(1369)
                    );
                }
            }
        }
    }

    #[test]
    fn expanded_closed_scalar_powers_keep_the_multinomial_size() {
        test_support::test_initialize();
        let dot = closed_dot(
            99817,
            [
                spenso::vector_symbol!("scope_binomial_p"),
                spenso::vector_symbol!("scope_binomial_q"),
            ],
        );
        let x = Atom::var(symbol!("scope_binomial_x"));
        let source = SymbolicTensor::infer((&x + &dot).pow(12)).unwrap();
        for (var, via_poly) in [
            (None, false),
            (None, true),
            (Some(x.as_view()), false),
            (Some(x.as_view()), true),
        ] {
            let expanded = source.expanded(var, via_poly).unwrap();
            assert_eq!(expanded.structure, source.structure);
            // Unrolling the power into copies with distinct dummies would emit
            // 2^12 = 4096 terms which only canonization merges again.
            assert_eq!(summands(&expanded), 13);
            let admitted = SymbolicTensor::infer(expanded.expression.clone()).unwrap();
            assert_eq!(admitted.structure, source.structure);
            assert_eq!(summands(&expanded.canonize().unwrap()), 13);
        }
    }

    #[test]
    fn selected_expansion_preserves_unrelated_closed_coefficients() {
        test_support::test_initialize();
        let dot = closed_dot(
            99811,
            [
                spenso::vector_symbol!("scope_selection_p"),
                spenso::vector_symbol!("scope_selection_q"),
            ],
        );
        let x = Atom::var(symbol!("scope_selection_x"));
        let y = Atom::var(symbol!("scope_selection_y"));
        let coefficient = &y + dot;
        for exponent in [Atom::num(5), Atom::num(-3), Atom::num(3) / Atom::num(2)] {
            for expression in [coefficient.pow(&exponent), &x * coefficient.pow(&exponent)] {
                let source = SymbolicTensor::infer(expression).unwrap();
                for via_poly in [false, true] {
                    let result = source.expanded(Some(x.as_view()), via_poly).unwrap();
                    assert_unchanged(&result, &source);
                }
            }
        }
    }

    #[test]
    fn expansion_preserves_powers_beyond_the_backend_materialization_range() {
        test_support::test_initialize();
        let dot = closed_dot(
            99813,
            [
                spenso::vector_symbol!("scope_bound_p"),
                spenso::vector_symbol!("scope_bound_q"),
            ],
        );
        let x = Atom::var(symbol!("scope_bound_x"));
        let coefficient = &x + dot;
        // The backend leaves these powers alone, so no copy may be prepared.
        for exponent in [4_294_967_296_i64, -4_294_967_296_i64] {
            let source = SymbolicTensor::infer(coefficient.pow(Atom::num(exponent))).unwrap();
            for var in [None, Some(x.as_view())] {
                for via_poly in [false, true] {
                    assert_unchanged(&source.expanded(var, via_poly).unwrap(), &source);
                }
            }
        }
    }

    #[test]
    fn scope_selected_function_occurrences_still_expand_after_renaming() {
        test_support::test_initialize();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let [mu, nu] = [99893, 99895].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let heads = [
            spenso::vector_symbol!("scope_function_p"),
            spenso::vector_symbol!("scope_function_q"),
            spenso::vector_symbol!("scope_function_r"),
        ];
        let vector =
            |head: Symbol, index: &Atom| FunctionBuilder::new(head).add_arg(index).finish();
        let [p, q, r] = heads.map(|head| vector(head, &mu));
        let x = Atom::var(symbol!("scope_function_x"));
        let canonical = |value: &Atom| {
            SymbolicTensor::infer(value.clone())
                .unwrap()
                .canonize()
                .unwrap()
                .expression
        };
        // The selected `q(mu)` closes a pair inside the sum, whose name the
        // neighbouring `r(mu)` reuses. Separation renames the selected copy.
        let source = SymbolicTensor::infer((&x + &p * &q) * &r).unwrap();
        let expected = &x * &r + vector(heads[0], &nu) * vector(heads[1], &nu) * &r;
        for via_poly in [false, true] {
            let expanded = source.expanded(Some(q.as_view()), via_poly).unwrap();
            assert_eq!(summands(&expanded), 2);
            assert_eq!(expanded.structure, source.structure);
            assert_eq!(canonical(&expanded.expression), canonical(&expected));
        }
    }

    #[test]
    fn function_selectors_without_a_separated_scope_expand_as_the_backend() {
        test_support::test_initialize();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let [mu, nu] = [99897, 99899].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let vector =
            |head: Symbol, index: &Atom| FunctionBuilder::new(head).add_arg(index).finish();
        let p = vector(spenso::vector_symbol!("scope_backend_p"), &mu);
        let q = vector(spenso::vector_symbol!("scope_backend_q"), &mu);
        let r = |index| vector(spenso::vector_symbol!("scope_backend_r"), index);
        let s = vector(spenso::vector_symbol!("scope_backend_s"), &mu);
        let x = Atom::var(symbol!("scope_backend_x"));
        let y = Atom::var(symbol!("scope_backend_y"));
        let closed = &x + &p * &q;
        let cases = [
            // Selected inside the closed scope, with no colliding neighbour.
            (&closed * r(&nu), q.clone()),
            // Selected only outside: the open sum is distributed, and the
            // closed sum stays one factor of each term.
            (&closed * (&y * &s + r(&mu)), r(&mu)),
            // Selected only outside a colliding closed sum, which stays intact.
            (&closed * r(&mu), r(&mu)),
        ];
        for (case, (source, selector)) in cases.into_iter().enumerate() {
            let source = SymbolicTensor::infer(source).unwrap();
            for via_poly in [false, true] {
                let expanded = source.expanded(Some(selector.as_view()), via_poly).unwrap();
                let backend = if via_poly {
                    source
                        .expression
                        .expand_via_poly::<i32, Atom>(Some(selector.clone()))
                } else {
                    source.expression.expand_in(selector.as_view())
                };
                if case == 1 && via_poly {
                    assert!(SymbolicTensor::infer(backend).is_err());
                    let closed = &x
                        + vector(p.get_symbol().unwrap(), &nu)
                            * vector(q.get_symbol().unwrap(), &nu);
                    let expected = SymbolicTensor::infer(closed * (&y * &s + r(&mu)))
                        .unwrap()
                        .expanded(None, false)
                        .unwrap()
                        .canonize()
                        .unwrap();
                    assert_eq!(
                        expanded
                            .expanded(None, false)
                            .unwrap()
                            .canonize()
                            .unwrap()
                            .expression,
                        expected.expression,
                    );
                } else {
                    assert_eq!(expanded.expression, backend);
                }
                if case == 1 {
                    let terms = match expanded.expression.as_view() {
                        AtomView::Add(sum) => sum.iter().map(|term| term.to_owned()).collect(),
                        _ => vec![expanded.expression.clone()],
                    };
                    // Only the selected open sum is distributed: the closed sum
                    // stays one factor of each term. The polynomial route renames
                    // its pair, which would otherwise meet `s(mu)` and `r(mu)`.
                    assert_eq!(terms.len(), 2, "output={}", expanded.expression);
                    for term in &terms {
                        let AtomView::Mul(product) = term.as_view() else {
                            panic!("{term}");
                        };
                        let sums = product
                            .iter()
                            .filter(|factor| matches!(factor, AtomView::Add(_)))
                            .collect::<Vec<_>>();
                        assert_eq!(sums.len(), 1, "{term}");
                        if via_poly {
                            assert!(!sums[0].contains(mu.as_view()), "{term}");
                            assert_eq!(
                                SymbolicTensor::infer(sums[0].to_owned())
                                    .unwrap()
                                    .canonize()
                                    .unwrap()
                                    .expression,
                                SymbolicTensor::infer(closed.clone())
                                    .unwrap()
                                    .canonize()
                                    .unwrap()
                                    .expression,
                            );
                        } else {
                            assert_eq!(sums[0], closed.as_view());
                        }
                    }
                }
                assert_eq!(expanded.structure, source.structure);
                let admitted = SymbolicTensor::infer(expanded.expression.clone()).unwrap();
                assert_eq!(admitted.structure, source.structure);
            }
        }
    }

    #[test]
    fn function_selection_skips_occupied_and_conflicting_carrier_names() {
        test_support::test_initialize();
        // The payload contains a plain symbol spelled like the first carrier
        // candidate, and the second one exists with other attributes.
        let occupied = SymbolBuilder::new(NamespacedSymbol::parse("idenso::expansion_selector_0"))
            .build()
            .unwrap();
        let _ = SymbolBuilder::new(NamespacedSymbol::parse("idenso::expansion_selector_1"))
            .with_attributes(vec![SymbolAttribute::Symmetric])
            .build();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let [mu, nu] = [99901, 99903].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let heads = [
            spenso::vector_symbol!("scope_carrier_p"),
            spenso::vector_symbol!("scope_carrier_q"),
            spenso::vector_symbol!("scope_carrier_r"),
        ];
        let vector =
            |head: Symbol, index: &Atom| FunctionBuilder::new(head).add_arg(index).finish();
        let [p, q, r] = heads.map(|head| vector(head, &mu));
        let x = Atom::var(symbol!("scope_carrier_x"));
        let z = Atom::var(occupied);
        let canonical = |value: &Atom| {
            SymbolicTensor::infer(value.clone())
                .unwrap()
                .canonize()
                .unwrap()
                .expression
        };
        let source = SymbolicTensor::infer(&z * (&x + &p * &q) * &r).unwrap();
        // Marking with `z` would turn it into one when the marker is removed.
        let expected = &z * &x * &r + &z * vector(heads[0], &nu) * vector(heads[1], &nu) * &r;
        for via_poly in [false, true] {
            let expanded = source.expanded(Some(q.as_view()), via_poly).unwrap();
            assert_eq!(summands(&expanded), 2);
            assert_eq!(expanded.structure, source.structure);
            assert_eq!(canonical(&expanded.expression), canonical(&expected));
        }
    }

    #[test]
    fn function_selection_without_separation_replays_no_callbacks() {
        test_support::test_initialize();
        let payload_calls = Arc::new(AtomicUsize::new(0));
        let selector_calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&payload_calls);
        let hook = symbolica::symbol!(
            "scope_callback_hook",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let selected_seen = Arc::clone(&selector_calls);
        let selected_head = spenso::vector_symbol!(
            "scope_callback_q",
            norm = move |_, _| {
                selected_seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let [mu, nu] = [99905, 99907].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let vector =
            |head: Symbol, index: &Atom| FunctionBuilder::new(head).add_arg(index).finish();
        let p = vector(spenso::vector_symbol!("scope_callback_p"), &mu);
        let q = vector(selected_head, &mu);
        let r = vector(spenso::vector_symbol!("scope_callback_r"), &nu);
        let x = Atom::var(symbol!("scope_callback_x"));
        let y = Atom::var(symbol!("scope_callback_y"));
        let hooked = |value: &Atom| FunctionBuilder::new(hook).add_arg(value).finish();
        // The selected sum is distributed, but nothing collides with its pair.
        let source = SymbolicTensor::infer(hooked(&y) * (&x * hooked(&x) + &p * &q) * &r).unwrap();
        for via_poly in [false, true] {
            // The prepass is the only step which the PR head does not run. On
            // this input it observes without constructing, normalizing, or
            // allocating, so every callback runs exactly as at the PR head.
            payload_calls.store(0, Ordering::Relaxed);
            selector_calls.store(0, Ordering::Relaxed);
            let mut allocations = 0;
            let separated = ExplicitIndexOccurrences::separate_expansion_scopes(
                source.expression.as_view(),
                Some(q.as_view()),
                via_poly,
                &mut || {
                    allocations += 1;
                    AbstractIndex::Normal(0)
                },
            );
            assert!(separated.is_none());
            assert_eq!(allocations, 0);
            assert_eq!(payload_calls.load(Ordering::Relaxed), 0);
            assert_eq!(selector_calls.load(Ordering::Relaxed), 0);
            // The original payload and selector then reach the backend.
            let expanded = source.expanded(Some(q.as_view()), via_poly).unwrap();
            let backend = if via_poly {
                source
                    .expression
                    .expand_via_poly::<i32, Atom>(Some(q.clone()))
            } else {
                source.expression.expand_in(q.as_view())
            };
            assert_eq!(expanded.expression, backend);
        }
    }

    #[test]
    fn native_delta_normalization_is_intrinsic_for_fixed_components() {
        crate::test_support::test_initialize();
        assert!(InterfaceInference::intrinsic_normalization_head(ETS.delta));
        let mut inference = InterfaceInference::default();
        for dimension in [
            Atom::num(4),
            Atom::var(symbolica::symbol!("intrinsic_delta_d")),
        ] {
            for component in 0..4 {
                let selector = symbolica::function!(AIND_SYMBOLS.cind, component);
                let basis = symbolica::function!(
                    ETS.delta,
                    &selector,
                    spenso::mink!(
                        spenso::structure::dimension::Dimension::try_from(dimension.as_view())
                            .unwrap()
                    )
                );
                assert!(InterfaceInference::normalization_is_intrinsic(
                    basis.as_view()
                ));
                assert!(inference.normalization_is_intrinsic_cached(basis.as_view()));
                let same = symbolica::function!(ETS.delta, &selector, &selector);
                let different = symbolica::function!(
                    ETS.delta,
                    &selector,
                    symbolica::function!(AIND_SYMBOLS.cind, (component + 1) % 4)
                );
                assert_eq!(same, Atom::one());
                assert!(different.is_zero());
            }
        }
    }

    #[test]
    fn native_delta_normalization_keeps_descendant_callback_guards() {
        use spenso::structure::representation::BaseRepName;
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let hook = symbolica::symbol!("intrinsic_delta_payload"; Scalar;
            norm = move |_, _| { observed.fetch_add(1, Ordering::Relaxed); });
        let observed = Arc::clone(&calls);
        let custom = spenso::vector_symbol!(
            "intrinsic_delta_custom",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let payload = symbolica::function!(hook, 7);
        let selector = symbolica::function!(AIND_SYMBOLS.cind, 0);
        let nested = symbolica::function!(
            ETS.delta,
            &selector,
            symbolica::function!(
                spenso::vector_symbol!("intrinsic_delta_p"),
                &payload,
                spenso::mink!(4)
            )
        );
        let dimension = symbolica::function!(
            ETS.delta,
            &selector,
            symbolica::function!(Minkowski::selfless_symbol(), &payload)
        );
        let custom = symbolica::function!(custom, spenso::mink!(4));
        calls.store(0, Ordering::Relaxed);
        let mut inference = InterfaceInference::default();
        for value in [&nested, &dimension, &custom] {
            assert!(!InterfaceInference::normalization_is_intrinsic(
                value.as_view()
            ));
            assert!(!inference.normalization_is_intrinsic_cached(value.as_view()));
            assert!(!inference.rewrites_preserve_leaf_interfaces(value.as_view()));
        }
        assert_eq!(calls.load(Ordering::Relaxed), 0);
        let source = SymbolicTensor::infer(spenso::dot!(&custom, &custom).pow(2)).unwrap();
        calls.store(0, Ordering::Relaxed);
        assert!(
            source
                .coefficient_list(spenso::shadowing::TensorCollectFilter::<0>::Tensors)
                .is_err()
        );
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn structural_dual_normalization_is_intrinsic_without_hiding_nested_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        crate::test_support::test_initialize();
        let port = ColorAntiFundamental {}
            .new_rep(3)
            .slot::<AbstractIndex, _>(97111)
            .to_atom();
        assert!(
            matches!(port.as_view(), AtomView::Fun(function) if function.get_symbol() == AIND_SYMBOLS.dind)
        );
        let leaf = FunctionBuilder::new(spenso::tensor_symbol!("intrinsic_dual_leaf"))
            .add_arg(&port)
            .finish();
        assert!(InterfaceInference::intrinsic_normalization_head(
            AIND_SYMBOLS.dind
        ));
        assert!(InterfaceInference::normalization_is_intrinsic(
            leaf.as_view()
        ));
        let mut inference = InterfaceInference::default();
        assert!(inference.normalization_is_intrinsic_cached(leaf.as_view()));
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let hook = symbolica::symbol!(
            "intrinsic_dual_payload_hook",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let nested = FunctionBuilder::new(AIND_SYMBOLS.dind)
            .add_arg(FunctionBuilder::new(hook).add_arg(&port).finish())
            .finish();
        calls.store(0, Ordering::Relaxed);
        assert!(!InterfaceInference::normalization_is_intrinsic(
            nested.as_view()
        ));
        assert!(!inference.normalization_is_intrinsic_cached(nested.as_view()));
        assert_eq!(
            calls.load(Ordering::Relaxed),
            0,
            "observing the dual wrapper must neither skip nor execute its user payload"
        );
    }

    #[test]
    fn admitted_index_callbacks_are_opaque_to_interface_proofs_and_contraction() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let label = symbolica::symbol!(
            "intrinsic_index_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            },
            tags = [SPENSO_TAG.index.clone()]
        );
        let a = spenso::mink!(4, symbolica::function!(label, 7, 0));
        let b = spenso::mink!(4, symbolica::function!(label, 7, 1));
        let spectator = symbolica::parse_lit!((intrinsic_index_x + intrinsic_index_y) ^ 9);
        let vector = spenso::vector_symbol!("intrinsic_index_vector");
        let expected = &spectator * symbolica::function!(vector, &b);
        let source = SymbolicTensor::infer(
            &spectator * spenso::g!(&a, &b) * symbolica::function!(vector, &a),
        )
        .unwrap();
        calls.store(0, Ordering::Relaxed);
        let mut inference = InterfaceInference::default();
        assert!(InterfaceInference::normalization_is_intrinsic(
            source.expression.as_view()
        ));
        assert!(inference.normalization_is_intrinsic_cached(source.expression.as_view()));
        assert!(inference.rewrites_preserve_leaf_interfaces(source.expression.as_view()));
        assert_eq!(
            calls.load(Ordering::Relaxed),
            0,
            "interface proof replayed label callbacks"
        );
        let result = source.contract(Default::default()).unwrap();
        assert_eq!(result.expression, expected);
        assert_eq!(result.structure, source.structure);
        assert_eq!(
            result.reduction_status(),
            crate::tensor::ReductionStatus::Complete
        );
        let again = result.contract(Default::default()).unwrap();
        assert_eq!(again.expression, expected);
        assert_eq!(again.structure, result.structure);
        assert_eq!(
            again.reduction_status(),
            crate::tensor::ReductionStatus::Complete
        );
        assert_eq!(
            calls.load(Ordering::Relaxed),
            0,
            "copying index identities must not invoke their normalizers"
        );

        // The same function outside an index is ordinary metadata and remains
        // a callback boundary. Its syntax alone is not a general exemption.
        let metadata = symbolica::function!(label, 9, 0);
        let with_metadata = symbolica::function!(vector, metadata, &b);
        calls.store(0, Ordering::Relaxed);
        assert!(!InterfaceInference::normalization_is_intrinsic(
            with_metadata.as_view()
        ));
        assert!(!inference.normalization_is_intrinsic_cached(with_metadata.as_view()));
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn structural_cyclic_normalization_checks_nested_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        crate::test_support::test_initialize();
        let a = Atom::var(symbolica::symbol!("intrinsic_cyclic_a"));
        let b = Atom::var(symbolica::symbol!("intrinsic_cyclic_b"));
        let cyclic = shadowing::cyclic([&a, &b]);
        assert!(matches!(cyclic.as_view(), AtomView::Fun(function)
            if function.get_symbol() == *shadowing::CYCLIC));
        assert!(InterfaceInference::intrinsic_normalization_head(
            *shadowing::CYCLIC
        ));
        let mut inference = InterfaceInference::default();
        assert!(InterfaceInference::normalization_is_intrinsic(
            cyclic.as_view()
        ));
        assert!(inference.normalization_is_intrinsic_cached(cyclic.as_view()));
        let symmetric = shadowing::sym([&a, &b]);
        assert_eq!(shadowing::cyclic([&symmetric]), symmetric);

        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let hook = symbolica::symbol!(
            "intrinsic_cyclic_payload_hook",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        assert!(!InterfaceInference::intrinsic_normalization_head(hook));
        let payload = FunctionBuilder::new(hook).add_arg(&a).finish();
        let nested = shadowing::cyclic([&payload, &b]);
        let nested_symmetric = shadowing::sym([&payload, &b]);
        let unwrapped = shadowing::cyclic([&nested_symmetric]);
        assert_eq!(unwrapped, nested_symmetric);
        calls.store(0, Ordering::Relaxed);
        for expression in [&nested, &unwrapped] {
            assert!(!InterfaceInference::normalization_is_intrinsic(
                expression.as_view()
            ));
            assert!(!inference.normalization_is_intrinsic_cached(expression.as_view()));
            assert!(!inference.rewrites_preserve_leaf_interfaces(expression.as_view()));
        }
        assert_eq!(
            calls.load(Ordering::Relaxed),
            0,
            "observing cyclic and unwrapped sym must neither skip nor execute user payloads"
        );
    }

    #[test]
    fn intrinsic_subtree_cache_reuses_complete_normalized_functions() {
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let p = FunctionBuilder::new(spenso::vector_symbol!("intrinsic_subtree_p"))
            .add_arg(rep.to_symbolic([]))
            .finish();
        let metadata = FunctionBuilder::new(symbolica::symbol!("intrinsic_subtree_meta"; Scalar))
            .add_arg(&p)
            .finish();
        let leaf = FunctionBuilder::new(spenso::tensor_symbol!("intrinsic_subtree_tensor"))
            .add_arg(metadata)
            .add_arg(&p)
            .add_arg(
                rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(731))
                    .to_atom(),
            )
            .finish();
        assert!(leaf.as_view().get_byte_size() <= InterfaceInference::CACHE_KEY_BYTES);
        let repeated = Atom::var(symbolica::symbol!("intrinsic_subtree_x")) * &leaf
            + Atom::var(symbolica::symbol!("intrinsic_subtree_y")) * &leaf;
        let mut inference = InterfaceInference::default();
        INTRINSIC_HEAD_CHECKS.with(|n| n.set(0));
        assert!(InterfaceInference::normalization_is_intrinsic(
            leaf.as_view()
        ));
        let once = INTRINSIC_HEAD_CHECKS.with(|n| n.get());
        assert!(once > 3);
        INTRINSIC_HEAD_CHECKS.with(|n| n.set(0));
        LEAF_PROOF_CALLS.with(|n| n.set(0));
        assert!(inference.rewrites_preserve_leaf_interfaces(repeated.as_view()));
        assert_eq!(INTRINSIC_HEAD_CHECKS.with(|n| n.get()), once);
        assert_eq!(LEAF_PROOF_CALLS.with(|n| n.get()), 1);
        assert_eq!(inference.intrinsic_functions.len(), 1);
        INTRINSIC_HEAD_CHECKS.with(|n| n.set(0));
        assert!(inference.rewrites_preserve_leaf_interfaces(repeated.as_view()));
        assert!(inference.algebra_preserves_leaf_interfaces(leaf.as_view()));
        assert_eq!(INTRINSIC_HEAD_CHECKS.with(|n| n.get()), 0);
        assert_eq!(LEAF_PROOF_CALLS.with(|n| n.get()), 1);
    }

    #[test]
    fn intrinsic_subtree_cache_preserves_global_barrier_and_failure_tail() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let m = spenso::mink!(4);
        let compact = FunctionBuilder::new(spenso::vector_symbol!("intrinsic_subtree_barrier_p"))
            .add_arg(&m)
            .finish();
        let head = spenso::tensor_symbol!("intrinsic_subtree_barrier_tensor");
        let plain_metadata =
            FunctionBuilder::new(symbolica::symbol!("intrinsic_subtree_barrier_meta"; Scalar))
                .add_arg(&compact)
                .finish();
        let leaf = FunctionBuilder::new(head)
            .add_arg(plain_metadata)
            .add_arg(spenso::mink!(4, 739))
            .finish();
        let seen = Arc::clone(&calls);
        let hook = symbolica::symbol!(
            "intrinsic_subtree_hook",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let bad = FunctionBuilder::new(hook).add_arg(&compact).finish();
        let tail = FunctionBuilder::new(symbolica::symbol!("intrinsic_subtree_tail"))
            .add_arg(&compact)
            .finish();
        // Products order distinct function heads by Symbol ID. Verify that
        // the safe structured leaf precedes rejection and a sibling follows.
        let source = &leaf * &bad * &tail;
        let AtomView::Mul(terms) = source.as_view() else {
            unreachable!()
        };
        assert_eq!(
            terms.iter().collect::<Vec<_>>(),
            vec![leaf.as_view(), bad.as_view(), tail.as_view()]
        );
        calls.store(0, Ordering::Relaxed);
        INTRINSIC_HEAD_CHECKS.with(|n| n.set(0));
        assert!(!InterfaceInference::normalization_is_intrinsic(
            source.as_view()
        ));
        let old_checks = INTRINSIC_HEAD_CHECKS.with(|n| n.get());
        let mut inference = InterfaceInference::default();
        LEAF_PROOF_CALLS.with(|n| n.set(0));
        INTRINSIC_HEAD_CHECKS.with(|n| n.set(0));
        assert!(!inference.rewrites_preserve_leaf_interfaces(source.as_view()));
        assert_eq!(INTRINSIC_HEAD_CHECKS.with(|n| n.get()), old_checks);
        assert_eq!(LEAF_PROOF_CALLS.with(|n| n.get()), 0);
        // A prior positive algebra proof must not bypass the global barrier.
        assert!(inference.algebra_preserves_leaf_interfaces(leaf.as_view()));
        LEAF_PROOF_CALLS.with(|n| n.set(0));
        for _ in 0..2 {
            assert!(!inference.rewrites_preserve_leaf_interfaces(source.as_view()));
            assert!(
                !inference
                    .intrinsic_functions
                    .contains(bad.as_view().get_data())
            );
        }
        assert_eq!(LEAF_PROOF_CALLS.with(|n| n.get()), 0);
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn intrinsic_subtree_cache_keeps_raw_negative_and_capacity_boundaries() {
        use symbolica::atom::UserData;
        let head = symbolica::symbol!("intrinsic_subtree_capacity");
        let mut inference = InterfaceInference::default();
        let mut first = None;
        for i in 0..InterfaceInference::CACHE_ENTRIES + 2 {
            let leaf = FunctionBuilder::new(head)
                .add_arg(Atom::num(i as i64))
                .finish();
            assert!(inference.normalization_is_intrinsic_cached(leaf.as_view()));
            first.get_or_insert_with(|| leaf.clone());
        }
        assert_eq!(
            inference.intrinsic_functions.len(),
            InterfaceInference::CACHE_ENTRIES
        );
        let first = first.unwrap();
        INTRINSIC_HEAD_CHECKS.with(|n| n.set(0));
        assert!(inference.normalization_is_intrinsic_cached(first.as_view()));
        assert_eq!(INTRINSIC_HEAD_CHECKS.with(|n| n.get()), 0);
        let missing = FunctionBuilder::new(head).add_arg(Atom::num(1000)).finish();
        INTRINSIC_HEAD_CHECKS.with(|n| n.set(0));
        for _ in 0..2 {
            assert!(inference.normalization_is_intrinsic_cached(missing.as_view()));
        }
        assert_eq!(INTRINSIC_HEAD_CHECKS.with(|n| n.get()), 2);
        // A raw view of a previously admitted call must still traverse its head.
        let mut encoded = Vec::new();
        UserData::Atom(first.clone()).write(&mut encoded).unwrap();
        let mut raw_head = Atom::new();
        raw_head.to_fun(head);
        let offset = encoded.len() - first.as_view().get_byte_size();
        encoded[offset] = raw_head.as_view().get_data()[0];
        let UserData::Atom(raw) = UserData::read(&mut encoded.as_slice()).unwrap() else {
            unreachable!()
        };
        assert!(raw.as_view().needs_normalization());
        INTRINSIC_HEAD_CHECKS.with(|n| n.set(0));
        for _ in 0..2 {
            assert!(inference.normalization_is_intrinsic_cached(raw.as_view()));
        }
        assert_eq!(INTRINSIC_HEAD_CHECKS.with(|n| n.get()), 2);
        assert!(
            !inference
                .intrinsic_functions
                .contains(raw.as_view().get_data())
        );
        let oversized = FunctionBuilder::new(head)
            .add_args((0..300).map(Atom::num))
            .finish();
        assert!(oversized.as_view().get_byte_size() > InterfaceInference::CACHE_KEY_BYTES);
        let mut fresh = InterfaceInference::default();
        for _ in 0..2 {
            assert!(fresh.normalization_is_intrinsic_cached(oversized.as_view()));
        }
        assert!(fresh.intrinsic_functions.is_empty());
    }

    #[test]
    fn intrinsic_subtree_cache_matches_two_phase_proof_for_reused_owners() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        use symbolica::atom::{EvaluationInfo, SymbolBuilder};
        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let hook = symbolica::symbol!(
            "intrinsic_subtree_matrix_hook",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let seen = Arc::clone(&calls);
        let eval = SymbolBuilder::new(symbolica::wrap_symbol!("intrinsic_subtree_matrix_eval"))
            .with_evaluation_info(EvaluationInfo::new().register(move |args: &[f64]| {
                seen.fetch_add(1, Ordering::Relaxed);
                args[0]
            }))
            .build()
            .unwrap();
        let tensor = spenso::tensor_symbol!("intrinsic_subtree_matrix_tensor");
        let m = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let compact = FunctionBuilder::new(spenso::vector_symbol!("intrinsic_subtree_matrix_p"))
            .add_arg(m.to_symbolic([]))
            .finish();
        let slot = m
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(745))
            .to_atom();
        let leaf = FunctionBuilder::new(tensor).add_arg(&slot).finish();
        let dot = FunctionBuilder::new(ETS.metric)
            .add_arg(&compact)
            .add_arg(&compact)
            .finish();
        let mut values = vec![
            Atom::num(0),
            leaf.clone(),
            dot.clone(),
            dot.pow(Atom::num(3)),
            leaf.pow(Atom::num(3)),
            FunctionBuilder::new(tensor)
                .add_arg(&compact)
                .add_arg(&slot)
                .finish(),
            FunctionBuilder::new(tensor)
                .add_arg(m.to_symbolic([]))
                .finish(),
            FunctionBuilder::new(tensor)
                .add_arg(
                    m.slot::<AbstractIndex, _>(AbstractIndex::Open {
                        owner: 746,
                        axis: 0,
                    })
                    .to_atom(),
                )
                .finish(),
            FunctionBuilder::new(ETS.metric).add_arg(&slot).finish(),
            FunctionBuilder::new(
                spenso::tensor_symbol!("intrinsic_subtree_matrix_symmetric"; Symmetric),
            )
            .add_arg(&slot)
            .finish(),
        ];
        for head in [hook, eval] {
            let bad = FunctionBuilder::new(head).add_arg(Atom::num(2)).finish();
            let opaque =
                FunctionBuilder::new(symbolica::symbol!("intrinsic_subtree_matrix_scalar"; Scalar))
                    .add_arg(&bad)
                    .finish();
            values.push(
                FunctionBuilder::new(tensor)
                    .add_arg(&opaque)
                    .add_arg(&slot)
                    .finish(),
            );
            values.push(opaque);
            values.push(Atom::var(symbolica::symbol!("intrinsic_subtree_matrix_x")).pow(&bad));
        }
        calls.store(0, Ordering::Relaxed);
        for policy in [
            LeafInference::Materialize,
            LeafInference::Observe,
            LeafInference::ObserveEncoded,
        ] {
            let mut inference = InterfaceInference {
                leaf_inference: policy,
                ..Default::default()
            };
            for value in values.iter().chain(values.iter().rev()) {
                let bytes = value.as_view().get_data().to_vec();
                let expected = InterfaceInference::normalization_is_intrinsic(value.as_view())
                    && InterfaceInference {
                        leaf_inference: policy,
                        ..Default::default()
                    }
                    .algebra_preserves_leaf_interfaces(value.as_view());
                assert_eq!(
                    inference.rewrites_preserve_leaf_interfaces(value.as_view()),
                    expected,
                    "{value}"
                );
                assert_eq!(value.as_view().get_data(), bytes);
                assert_eq!(calls.load(Ordering::Relaxed), 0);
            }
        }
    }

    fn explicit_interface(
        representation: Representation<LibraryRep>,
        index: AbstractIndex,
    ) -> PartialStructure {
        PartialStructure::from_logical_slots([representation.slot(PartialIndex::Explicit(index))])
    }

    #[test]
    fn checked_results_validate_scopes_once_and_preserve_logical_ports() {
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let head = spenso::tensor_symbol!("checked_result_layout");
        for slots in [
            [
                euclidean.slot(PartialIndex::Explicit(AbstractIndex::Normal(74001))),
                minkowski.slot(PartialIndex::Explicit(AbstractIndex::Normal(74003))),
            ],
            [
                euclidean.slot(PartialIndex::open(0)),
                minkowski.slot(PartialIndex::open(1)),
            ],
        ] {
            let leaf = FunctionBuilder::new(head)
                .add_args(slots.map(composition::port_atom))
                .finish();
            let logical = if matches!(slots[0].aind, PartialIndex::Explicit(_)) {
                // Explicit ports can have a public order distinct from syntax.
                PartialStructure::from_logical_slots(slots.into_iter().rev())
            } else {
                PartialStructure::from_logical_slots(slots)
            };
            let source = SymbolicTensor::<PartialStructure>::new(leaf.clone(), logical).unwrap();
            let expression = &leaf * Atom::num(2);
            SCOPE_VALIDATIONS.with(|count| count.set(0));
            SymbolicTensor::validate_interface(&expression, &source.structure).unwrap();
            let reference =
                SymbolicTensor::from_produced_parts(expression.clone(), source.structure.clone())
                    .unwrap();
            // Only the explicit validation above admits raw input.
            assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 1);

            SCOPE_VALIDATIONS.with(|count| count.set(0));
            let result = source.with_checked_expression(expression).unwrap();
            assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 1);
            assert_eq!(result, reference);
            assert_eq!(
                result.structure.logical_slots(),
                source.structure.logical_slots()
            );

            let zero = source.with_checked_expression(Atom::Zero).unwrap();
            assert!(zero.expression.is_zero());
            assert_eq!(zero.structure, source.structure);
        }
    }

    #[test]
    fn checked_results_preserve_error_precedence_and_branch_validation() {
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let a = representation.slot::<AbstractIndex, _>(AbstractIndex::Normal(74011));
        let b = representation.slot::<AbstractIndex, _>(AbstractIndex::Normal(74013));
        let head = spenso::tensor_symbol!("checked_result_errors");
        let leaf = |slots: &[Atom]| FunctionBuilder::new(head).add_args(slots).finish();
        let tensor = SymbolicTensor::<PartialStructure>::new(
            leaf(&[a.to_atom()]),
            explicit_interface(representation, a.aind()),
        )
        .unwrap();
        let scalar = SymbolicTensor::<PartialStructure>::new(
            Atom::one(),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        let placeholder = Atom::var(SPENSO_TAG.chain_in);
        for (source, expression, expected) in [
            (&tensor, placeholder.clone(), "compatible tensor interface"),
            (&scalar, placeholder, "chain"),
            (
                &tensor,
                leaf(&[a.to_atom()]) + leaf(&[b.to_atom()]),
                "compatible",
            ),
            (
                &scalar,
                leaf(&[a.to_atom(), a.to_atom(), a.to_atom()]),
                "more than two",
            ),
        ] {
            let reference = SymbolicTensor::validate_interface(&expression, &source.structure)
                .and_then(|()| {
                    SymbolicTensor::from_produced_parts(
                        expression.clone(),
                        source.structure.clone(),
                    )
                })
                .unwrap_err()
                .to_string();
            let actual = source
                .with_checked_expression(expression)
                .unwrap_err()
                .to_string();
            assert_eq!(actual, reference);
            assert!(actual.contains(expected), "{actual}");
        }
    }

    #[test]
    fn observed_network_indices_preserve_partial_observation_policy() {
        let m = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let a = m.slot::<AbstractIndex, _>(AbstractIndex::Open {
            owner: 101,
            axis: 0,
        });
        let b = m.slot::<AbstractIndex, _>(AbstractIndex::Open {
            owner: 102,
            axis: 0,
        });
        let tensor = spenso::tensor_symbol!("encoded_observation_tensor");
        for (atom, slots) in [
            (
                FunctionBuilder::new(tensor)
                    .add_arg(a.to_atom())
                    .add_arg(b.to_atom())
                    .finish(),
                vec![a, b],
            ),
            (
                FunctionBuilder::new(ETS.metric)
                    .add_arg(a.to_atom())
                    .add_arg(b.to_atom())
                    .finish(),
                vec![a, b],
            ),
        ] {
            let encoded = PartialStructure::from_logical_slots(
                slots
                    .iter()
                    .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
            );
            SymbolicTensor::validate_encoded_interface(&atom, &encoded).unwrap();
            let mut changed = slots.clone();
            changed[0].set_aind(AbstractIndex::Open {
                owner: 103,
                axis: 0,
            });
            let changed = PartialStructure::from_logical_slots(
                changed
                    .iter()
                    .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
            );
            assert!(SymbolicTensor::validate_encoded_interface(&atom, &changed).is_err());
            let unresolved = PartialStructure::from_logical_slots(
                slots
                    .iter()
                    .map(|slot| slot.rep().slot(PartialIndex::open(0))),
            );
            SymbolicTensor::validate_interface(&atom, &unresolved).unwrap();
        }
        // A bare slot is scalar metadata, not a standalone tensor leaf.
        let scalar = PartialStructure::from_logical_slots([]);
        SymbolicTensor::validate_encoded_interface(&a.to_atom(), &scalar).unwrap();
        assert!(
            SymbolicTensor::validate_encoded_interface(
                &a.to_atom(),
                &explicit_interface(m, a.aind())
            )
            .is_err()
        );
    }

    #[test]
    fn compact_bound_ports_preserve_surviving_logical_order() {
        let tensor = spenso::tensor_symbol!("bound_port_order_tensor");
        let vector = spenso::vector_symbol!("bound_port_order_vector");
        let m = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let e = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let slots = [
            m.slot(AbstractIndex::Normal(1)),
            e.slot(AbstractIndex::Normal(2)),
            m.slot(AbstractIndex::Normal(3)),
        ];
        for bound in 0..3 {
            let compact = FunctionBuilder::new(vector)
                .add_arg(slots[bound].rep().to_symbolic([]))
                .finish();
            let atom = FunctionBuilder::new(tensor)
                .add_args(slots.iter().enumerate().map(|(position, slot)| {
                    if position == bound {
                        compact.clone()
                    } else {
                        slot.to_atom()
                    }
                }))
                .finish();
            let expected = PartialStructure::from_logical_slots(
                slots
                    .iter()
                    .enumerate()
                    .filter(|(position, _)| *position != bound)
                    .map(|(_, slot)| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
            );
            let inferred = InterfaceInference::default()
                .infer_validated(atom.as_view())
                .unwrap();
            assert_eq!(inferred, expected, "bound argument {bound}");
            SymbolicTensor::<PartialStructure>::validate_interface(&atom, &expected).unwrap();
        }
        let compact = FunctionBuilder::new(vector)
            .add_arg(m.to_symbolic([]))
            .finish();
        let open = FunctionBuilder::new(tensor)
            .add_arg(e.to_symbolic([]))
            .add_arg(&compact)
            .add_arg(m.to_symbolic([]))
            .finish();
        let inferred = InterfaceInference::default()
            .infer_validated(open.as_view())
            .unwrap();
        assert_eq!(
            inferred
                .logical_slots()
                .iter()
                .map(IsAbstractSlot::rep)
                .collect::<Vec<_>>(),
            vec![e, m]
        );
        assert_eq!(inferred.open_positions(), vec![0, 1]);
        let repeated = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(&open)
            .add_arg(&open)
            .finish();
        let repeated = InterfaceInference::default()
            .infer_validated(repeated.as_view())
            .unwrap();
        assert_eq!(repeated.open_positions(), vec![0, 1, 2, 3]);
        let scalar = FunctionBuilder::new(tensor)
            .add_arg(&compact)
            .add_arg(&compact)
            .finish();
        assert!(
            InterfaceInference::default()
                .infer_validated(scalar.as_view())
                .unwrap()
                .canonical()
                .is_scalar()
        );
    }

    #[test]
    fn compact_bound_ports_check_metadata_without_replaying_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let outer_calls = Arc::clone(&calls);
        let tensor = spenso::tensor_symbol!(
            "bound_port_callback_tensor",
            norm = move |_, _| {
                outer_calls.fetch_add(1, Ordering::Relaxed);
            }
        );
        let inner_calls = Arc::clone(&calls);
        let vector = spenso::vector_symbol!(
            "bound_port_callback_vector",
            norm = move |_, _| {
                inner_calls.fetch_add(1, Ordering::Relaxed);
            }
        );
        let m = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let metadata =
            FunctionBuilder::new(symbolica::symbol!("bound_port_scalar_metadata"; Scalar))
                .add_arg(m.to_symbolic([]))
                .finish();
        let compact = FunctionBuilder::new(vector)
            .add_arg(&metadata)
            .add_arg(m.to_symbolic([]))
            .finish();
        let atom = FunctionBuilder::new(tensor)
            .add_arg(&metadata)
            .add_arg(
                m.slot::<AbstractIndex, _>(AbstractIndex::Normal(17))
                    .to_atom(),
            )
            .add_arg(&compact)
            .add_arg(m.to_symbolic([]))
            .finish();
        calls.store(0, Ordering::Relaxed);
        let interface = InterfaceInference::default()
            .infer_validated(atom.as_view())
            .unwrap();
        assert_eq!(interface.canonical().order(), 2);
        assert_eq!(interface.open_positions(), vec![1]);
        SymbolicTensor::<PartialStructure>::validate_interface(&atom, &interface).unwrap();
        assert_eq!(calls.load(Ordering::Relaxed), 0);
        assert!(!InterfaceInference::default().algebra_preserves_leaf_interfaces(atom.as_view()));
        assert!(
            SymbolicTensor::<PartialStructure>::validate_interface(&Atom::one(), &interface)
                .is_err()
        );
        SymbolicTensor::<PartialStructure>::validate_interface(&Atom::Zero, &interface).unwrap();
        let scalar_atom = FunctionBuilder::new(tensor).add_arg(&compact).finish();
        let scalar = InterfaceInference::default()
            .infer_validated(scalar_atom.as_view())
            .unwrap();
        assert!(scalar.canonical().is_scalar());
    }

    #[test]
    fn compact_bound_ports_do_not_admit_tensor_metadata_or_noncompact_vectors() {
        let tensor = spenso::tensor_symbol!("bound_port_invalid_tensor");
        let vector = spenso::vector_symbol!("bound_port_invalid_vector");
        let m = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let rep = m.to_symbolic([]);
        let slot = m
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(3))
            .to_atom();
        let rank_two = FunctionBuilder::new(tensor)
            .add_arg(&rep)
            .add_arg(&rep)
            .finish();
        for invalid in [
            rank_two.clone(),
            FunctionBuilder::new(vector).add_arg(&slot).finish(),
            FunctionBuilder::new(vector)
                .add_arg(&rep)
                .add_arg(&rep)
                .finish(),
            FunctionBuilder::new(vector).add_arg(Atom::num(3)).finish(),
            FunctionBuilder::new(vector)
                .add_arg(&rank_two)
                .add_arg(&rep)
                .finish(),
        ] {
            let atom = FunctionBuilder::new(tensor)
                .add_arg(invalid)
                .add_arg(&slot)
                .finish();
            assert!(
                InterfaceInference::default()
                    .infer_validated(atom.as_view())
                    .is_err(),
                "{atom}"
            );
        }
        let scalar =
            FunctionBuilder::new(symbolica::symbol!("bound_port_invalid_scalar_wrapper"; Scalar))
                .add_arg(rank_two)
                .finish();
        let compact = FunctionBuilder::new(vector)
            .add_arg(&scalar)
            .add_arg(rep)
            .finish();
        let accepted = FunctionBuilder::new(tensor)
            .add_arg(scalar)
            .add_arg(slot)
            .add_arg(compact)
            .finish();
        assert_eq!(
            InterfaceInference::default()
                .infer_validated(accepted.as_view())
                .unwrap()
                .canonical()
                .order(),
            1
        );
    }

    #[test]
    fn intrinsic_plain_and_bound_leaf_results_reuse_the_validated_interface() {
        let m = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let tensor = spenso::tensor_symbol!("bound_reuse_tensor");
        let vector = spenso::vector_symbol!("bound_reuse_vector");
        let compact = FunctionBuilder::new(vector)
            .add_arg(m.to_symbolic([]))
            .finish();
        let momentum = FunctionBuilder::new(symbolica::symbol!("bound_reuse_momentum"; Scalar))
            .add_arg(&compact)
            .finish();
        let slot = m
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(7))
            .to_atom();
        let leaf = FunctionBuilder::new(tensor)
            .add_arg(&momentum)
            .add_arg(&compact)
            .add_arg(&slot)
            .finish();
        let unbound = FunctionBuilder::new(tensor)
            .add_arg(&momentum)
            .add_arg(&slot)
            .finish();
        let x = Atom::var(symbolica::symbol!("bound_reuse_x"));
        let y = Atom::var(symbolica::symbol!("bound_reuse_y"));
        for leaf in [leaf, unbound] {
            let original = SymbolicTensor::<PartialStructure>::infer(leaf.clone()).unwrap();
            assert!(
                InterfaceInference::default().rewrites_preserve_leaf_interfaces(leaf.as_view())
            );
            let rewritten = Atom::num(2) * &leaf;
            INFERENCE_CALLS.with(|count| count.set(0));
            INFERENCE_VISITS.with(|count| count.set(0));
            let result = original
                .with_rewritten_expression(rewritten.clone())
                .unwrap();
            assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
            assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
            assert_eq!(result.expression, rewritten);
            assert_eq!(result.structure, original.structure);
            let factored = SymbolicTensor::<PartialStructure>::infer((&x + &y) * &leaf).unwrap();
            let expanded = factored.expression.expand();
            assert_ne!(expanded, factored.expression);
            INFERENCE_CALLS.with(|count| count.set(0));
            INFERENCE_VISITS.with(|count| count.set(0));
            let result = factored.with_algebra_result(expanded.clone()).unwrap();
            assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
            assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
            assert_eq!(result.expression, expanded);
            assert_eq!(result.structure, factored.structure);
        }
        let unresolved = FunctionBuilder::new(tensor)
            .add_arg(&compact)
            .add_arg(m.to_symbolic([]))
            .finish();
        assert!(
            !InterfaceInference::default().algebra_preserves_leaf_interfaces(unresolved.as_view())
        );
        let symmetric =
            FunctionBuilder::new(spenso::tensor_symbol!("bound_reuse_symmetric"; Symmetric))
                .add_arg(&compact)
                .add_arg(slot)
                .finish();
        assert!(
            !InterfaceInference::default().algebra_preserves_leaf_interfaces(symmetric.as_view())
        );
    }

    #[test]
    fn bound_leaf_proof_reuses_successes_but_rechecks_unsupported_leaves() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };

        let m = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let tensor = spenso::tensor_symbol!("bound_proof_cache_tensor");
        let vector = spenso::vector_symbol!("bound_proof_cache_vector");
        let compact = FunctionBuilder::new(vector)
            .add_arg(m.to_symbolic([]))
            .finish();
        let metadata =
            FunctionBuilder::new(symbolica::symbol!("bound_proof_cache_metadata"; Scalar))
                .add_arg(&compact)
                .finish();
        let slot = m
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(7))
            .to_atom();
        let leaf = FunctionBuilder::new(tensor)
            .add_arg(metadata)
            .add_arg(&compact)
            .add_arg(&slot)
            .finish();
        assert!(leaf.as_view().get_byte_size() <= InterfaceInference::CACHE_KEY_BYTES);
        let repeated = Atom::var(symbolica::symbol!("bound_proof_cache_x")) * &leaf
            + Atom::var(symbolica::symbol!("bound_proof_cache_y")) * &leaf;
        assert_eq!(repeated.nterms(), 2);
        let mut inference = InterfaceInference::default();
        LEAF_PROOF_CALLS.with(|count| count.set(0));
        assert!(inference.algebra_preserves_leaf_interfaces(repeated.as_view()));
        assert_eq!(LEAF_PROOF_CALLS.with(|count| count.get()), 1);
        assert!(inference.rewrites_preserve_leaf_interfaces(repeated.as_view()));
        assert_eq!(LEAF_PROOF_CALLS.with(|count| count.get()), 1);

        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let callback = spenso::tensor_symbol!(
            "bound_proof_cache_callback",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let callback = FunctionBuilder::new(callback)
            .add_arg(&compact)
            .add_arg(slot)
            .finish();
        let unresolved = FunctionBuilder::new(tensor)
            .add_arg(&compact)
            .add_arg(m.to_symbolic([]))
            .finish();
        calls.store(0, Ordering::Relaxed);
        LEAF_PROOF_CALLS.with(|count| count.set(0));
        for unsupported in [callback, unresolved] {
            for _ in 0..2 {
                assert!(!inference.algebra_preserves_leaf_interfaces(unsupported.as_view()));
            }
        }
        assert_eq!(LEAF_PROOF_CALLS.with(|count| count.get()), 4);
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn bound_leaf_reuse_rejects_callbacks_before_metadata_inference() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let metadata = symbolica::symbol!("bound_reuse_callback_metadata"; Scalar; norm=move|_,_|{seen.fetch_add(1,Ordering::Relaxed);});
        let seen = Arc::clone(&calls);
        let callback_vector = spenso::vector_symbol!(
            "bound_reuse_callback_operand",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let m = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let tensor = spenso::tensor_symbol!("bound_reuse_callback_container");
        let vector = spenso::vector_symbol!("bound_reuse_plain_operand");
        let compact = FunctionBuilder::new(vector)
            .add_arg(m.to_symbolic([]))
            .finish();
        let callback_compact = FunctionBuilder::new(callback_vector)
            .add_arg(m.to_symbolic([]))
            .finish();
        let meta = FunctionBuilder::new(metadata)
            .add_arg(Atom::num(3))
            .finish();
        let slot = m
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(9))
            .to_atom();
        for leaf in [
            FunctionBuilder::new(tensor)
                .add_arg(&meta)
                .add_arg(&slot)
                .finish(),
            FunctionBuilder::new(tensor)
                .add_arg(&meta)
                .add_arg(&compact)
                .add_arg(&slot)
                .finish(),
            FunctionBuilder::new(tensor)
                .add_arg(&callback_compact)
                .add_arg(&slot)
                .finish(),
        ] {
            calls.store(0, Ordering::Relaxed);
            assert!(
                !InterfaceInference::default().algebra_preserves_leaf_interfaces(leaf.as_view())
            );
            assert_eq!(calls.load(Ordering::Relaxed), 0);
        }
        crate::representations::initialize();
        let directed = ColorFundamental {}.new_rep(Dimension::Concrete(3));
        let directed_compact = FunctionBuilder::new(vector)
            .add_arg(directed.to_symbolic([]))
            .finish();
        let leaf = FunctionBuilder::new(tensor)
            .add_arg(directed_compact)
            .add_arg(slot)
            .finish();
        assert!(!InterfaceInference::default().algebra_preserves_leaf_interfaces(leaf.as_view()));
    }

    #[test]
    fn builtin_signature_errors_remain_structured_through_syntax_validation() {
        let atom = FunctionBuilder::new(AGS.gamma)
            .add_arg(Atom::one())
            .finish();
        assert!(matches!(
            InterfaceInference::default().infer_validated(atom.as_view()),
            Err(TensorInferenceError::InvalidBuiltinSignature {
                factory: "gamma",
                ..
            })
        ));
    }

    #[test]
    fn signature_construction_retains_logical_order_and_allocates_distinct_occurrences() {
        crate::representations::initialize();
        let signature = CS.t_strct(Dimension::Concrete(3), Dimension::Concrete(8));
        let first = SymbolicTensor::<PartialStructure>::from_signature(&signature).unwrap();
        let second = SymbolicTensor::<PartialStructure>::from_signature(&signature).unwrap();
        assert_ne!(first.expression, second.expression);
        assert_eq!(first.structure, second.structure);
        assert_eq!(first.rank(), 3);
        assert_eq!(first.structure.open_positions(), vec![0, 1, 2]);
        assert_eq!(
            first
                .structure
                .logical_slots()
                .iter()
                .map(IsAbstractSlot::rep)
                .collect::<Vec<_>>(),
            signature.layout().canonical_to_logical(
                &signature
                    .canonical()
                    .external_reps_iter()
                    .collect::<Vec<_>>()
            ),
        );

        let strips_ports = spenso::tensor_symbol!(
            "signature_constructor_strips_ports",
            norm = |_, output| {
                **output = Atom::one();
            }
        );
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let invalid = ExplicitKey::from_iter([representation], strips_ports, None);
        assert!(SymbolicTensor::<PartialStructure>::from_signature(&invalid).is_err());
    }

    #[test]
    fn signature_constructor_retains_raw_ports_for_network_lowering() {
        test_support::test_initialize();
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let head = tensor_symbol!("signature_retained_raw_ports");
        let signature = ExplicitKey::from_iter([euclidean, minkowski, minkowski], head, None);
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let value = SymbolicTensor::from_signature(&signature).unwrap();
        let AtomView::Fun(function) = value.expression.as_view() else {
            panic!("the generic signature keeps its constructor head");
        };
        let raw = value.proofs.leaf_interfaces.get().unwrap().lock().unwrap();
        assert_eq!(
            raw.get(value.expression.as_view().get_data())
                .unwrap()
                .1
                .iter()
                .map(Atom::as_view)
                .collect::<Vec<_>>(),
            function.iter().collect::<Vec<_>>()
        );
        drop(raw);
        let network = value
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(network.graph.dangling_indices().len(), 3);
        assert_eq!(value.structure.open_positions(), vec![0, 1, 2]);
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn materialized_leaf_uses_spenso_logical_order_through_scalar_metadata() {
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let slots = [
            rep.slot(AbstractIndex::from(31)),
            rep.slot(AbstractIndex::from(7)),
        ];
        let bundle = FunctionBuilder::new(spenso::structure::abstract_index::AIND_SYMBOLS.aind)
            .add_args(slots.map(|slot| slot.to_atom()))
            .finish();
        let head = spenso::tensor_symbol!("inference_bundle_layout");
        let invalid = FunctionBuilder::new(head).add_arg(bundle).finish();
        assert!(SymbolicTensor::<PartialStructure>::infer(invalid).is_err());
        // Scalar function metadata sends materializing inference through the
        // shared fallback without admitting a structural metadata bundle.
        let metadata = FunctionBuilder::new(SPENSO_TAG.scalar)
            .add_arg(Atom::num(3))
            .finish();
        assert!(matches!(metadata.as_view(), AtomView::Fun(_)));
        let atom = FunctionBuilder::new(head)
            .add_arg(metadata)
            .add_args(slots.map(|slot| slot.to_atom()))
            .finish();
        let value = SymbolicTensor::<PartialStructure>::infer(atom).unwrap();
        assert_eq!(
            value.structure.logical_slots(),
            slots.map(|slot| { slot.rep().slot(PartialIndex::Explicit(slot.aind())) })
        );
    }

    #[test]
    fn mixed_signature_layout_survives_physical_storage_order_and_a_noop_normalizer() {
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let plain = spenso::tensor_symbol!("mixed_signature_order");
        let callback = spenso::tensor_symbol!("mixed_signature_order_callback", norm = |_, _| {});
        for name in [plain, callback] {
            let signature = ExplicitKey::from_iter([euclidean, minkowski], name, None);
            let value = SymbolicTensor::<PartialStructure>::from_signature(&signature).unwrap();
            assert_eq!(
                value
                    .structure
                    .logical_slots()
                    .iter()
                    .map(IsAbstractSlot::rep)
                    .collect::<Vec<_>>(),
                vec![euclidean, minkowski],
            );
            assert_eq!(value.structure.open_positions(), vec![0, 1]);
            let dimension = Dimension::Concrete(6);
            let promoted = value.with_lorentz_dimension(dimension).unwrap();
            assert_eq!(
                promoted.structure.logical_slots(),
                vec![
                    euclidean.slot(PartialIndex::open(0)),
                    ExtendibleReps::MINKOWSKI
                        .new_rep(dimension)
                        .slot(PartialIndex::open(1)),
                ],
            );
            assert_eq!(
                promoted.expression,
                value
                    .expression
                    .with_lorentz_dimension(dimension.to_symbolic().as_view()),
            );
        }
    }

    #[test]
    fn lorentz_dimension_preserves_logical_ports_typed_zero_and_dispatch_flags() {
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let slots = [
            euclidean.slot(PartialIndex::open(0)),
            minkowski.slot(PartialIndex::Explicit(AbstractIndex::Normal(73301))),
            ExtendibleReps::MINKOWSKI
                .new_rep(Dimension::Concrete(3))
                .slot(PartialIndex::open(1)),
        ];
        let atom = FunctionBuilder::new(tensor_symbol!("dimension_logical_ports"))
            .add_arg(4)
            .add_args(slots.into_iter().rev().map(composition::port_atom))
            .finish();
        let dimension = Dimension::from(symbol!("shared_lorentz_dimension"));
        let expected_slots = [
            slots[0],
            ExtendibleReps::MINKOWSKI
                .new_rep(dimension)
                .slot(slots[1].aind),
            slots[2],
        ];
        for expression in [atom, Atom::Zero] {
            let value = SymbolicTensor::<PartialStructure>::from_produced_parts(
                expression,
                PartialStructure::from_logical_slots(slots),
            )
            .unwrap();
            INFERENCE_CALLS.with(|count| count.set(0));
            INFERENCE_VISITS.with(|count| count.set(0));
            let result = value.with_lorentz_dimension(dimension).unwrap();
            assert_eq!(result.structure.logical_slots(), expected_slots);
            assert_eq!(result.structure.open_positions(), vec![0, 2]);
            assert_eq!(result.is_metric, value.is_metric);
            assert_eq!(result.is_composite, value.is_composite);
            assert_eq!(result.expression.is_zero(), value.expression.is_zero());
            assert_eq!(result.with_lorentz_dimension(dimension).unwrap(), result);
            assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
            assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        }

        let gamma = SymbolicTensor::<PartialStructure>::from_signature(
            &AGS.gamma_strct::<AbstractIndex>(Dimension::Concrete(4)),
        )
        .unwrap();
        let promoted = gamma.with_lorentz_dimension(dimension).unwrap();
        assert_eq!(
            promoted
                .structure
                .logical_slots()
                .iter()
                .map(IsAbstractSlot::dim)
                .collect::<Vec<_>>(),
            vec![4.into(), 4.into(), dimension],
        );
        assert_eq!(
            promoted.structure.open_positions(),
            gamma.structure.open_positions()
        );

        let metric =
            SymbolicTensor::<PartialStructure>::from_produced_parts(
                FunctionBuilder::new(ETS.metric)
                    .add_args([73303, 73305].map(|index| {
                        minkowski
                            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                            .to_atom()
                    }))
                    .finish(),
                PartialStructure::from_logical_slots([73303, 73305].map(|index| {
                    minkowski.slot(PartialIndex::Explicit(AbstractIndex::Normal(index)))
                })),
            )
            .unwrap();
        let promoted = metric.with_lorentz_dimension(dimension).unwrap();
        assert!(promoted.is_metric);
        assert!(!promoted.is_composite);
    }

    #[test]
    fn identity_dimension_mapping_retains_reduction_state_and_zero_layout() {
        crate::test_support::test_initialize();
        let representation = ExtendibleReps::MINKOWSKI.new_rep(4);
        let ports = [73343, 73341]
            .map(|index| representation.slot(PartialIndex::Explicit(AbstractIndex::Normal(index))));
        let expression =
            FunctionBuilder::new(spenso::tensor_symbol!("dimension_noop_retained_matrix"))
                .add_args(ports.into_iter().rev().map(composition::port_atom))
                .finish();
        let value = SymbolicTensor::from_declared_parts(
            expression,
            PartialStructure::from_logical_slots(ports),
        )
        .unwrap();
        let d = Atom::var(symbolica::symbol!("dimension_noop_parameter"));
        let dimension = d.pow(2) - 1;
        let value = value
            .map_representation_dimensions(|_, _| dimension.clone())
            .unwrap();
        let zero = value.with_checked_expression(Atom::Zero).unwrap();
        for source in [value, zero] {
            let source = source.contract(Default::default()).unwrap();
            assert!(source.contraction_complete());
            let mut calls = Vec::new();
            let mapped = source
                .map_representation_dimensions(|representation, dimension| {
                    calls.push((representation, dimension.to_owned()));
                    dimension.to_owned()
                })
                .unwrap();
            assert_eq!(calls, vec![(representation.rep, dimension.clone())]);
            assert_eq!(mapped, source);
            assert!(mapped.same_reduction_state(&source));
            assert_eq!(
                mapped.structure.logical_slots(),
                source.structure.logical_slots()
            );
            assert!(
                mapped
                    .contract(Default::default())
                    .unwrap()
                    .same_reduction_state(&source)
            );

            // In particular, a zero's unchanged Atom cannot hide a dimension change.
            let specialized = source
                .map_representation_dimensions(|_, _| Atom::num(6))
                .unwrap();
            assert_eq!(specialized.is_zero(), source.is_zero());
            assert_eq!(
                specialized.structure.logical_slots(),
                ports.map(|slot| ExtendibleReps::MINKOWSKI.new_rep(6).slot(slot.aind))
            );
            assert!(!specialized.same_reduction_state(&source));
        }
    }

    #[test]
    fn identity_dimension_mapping_still_finishes_raw_parts() {
        test_support::test_initialize();
        let representation = ExtendibleReps::MINKOWSKI.new_rep(4);
        let slot = representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(73345))
            .to_atom();
        let malformed = FunctionBuilder::new(tensor_symbol!("dimension_noop_invalid_multiplicity"))
            .add_args([&slot, &slot, &slot])
            .finish();
        let scalar = PartialStructure::from_logical_slots([]);
        let expected = SymbolicTensor::from_produced_parts(malformed.clone(), scalar.clone())
            .unwrap_err()
            .to_string();
        let actual = SymbolicTensor::from_declared_parts(malformed, scalar).unwrap_err();
        assert_eq!(actual.to_string(), expected);

        // Interface admission can precede product finishing. Preserve that boundary.
        let x = Atom::var(symbol!("dimension_noop_unfinished_x"));
        let y = Atom::var(symbol!("dimension_noop_unfinished_y"));
        let expression = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(&x)
            .add_arg(&y)
            .finish();
        let admitted = SymbolicTensor::from_declared_parts(
            expression,
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        assert!(admitted.proofs.is_admitted());
        let finished = admitted
            .map_representation_dimensions(|_, dimension| dimension.to_owned())
            .unwrap();
        assert_eq!(finished.expression, x * y);
        assert!(finished.is_scalar());
    }

    #[test]
    fn mapped_color_dimensions_keep_encoded_ports_and_ranked_zero() {
        crate::test_support::test_initialize();
        let adjoint = ColorAdjoint {}.new_rep(8);
        let fundamental = ColorFundamental {}.new_rep(3);
        let metadata = symbolica::symbol!("mapped_color_scalar_metadata"; Scalar);
        let coefficient = FunctionBuilder::new(metadata)
            .add_arg(adjoint.to_symbolic([]))
            .finish();
        let expression = coefficient
            * FunctionBuilder::new(spenso::tensor_symbol!("mapped_color_dimensions"))
                .add_arg(fundamental.to_symbolic([]))
                .add_arg(adjoint.to_symbolic([Atom::var(symbolica::symbol!("mapped_color_index"))]))
                .add_arg(adjoint.to_symbolic([]))
                .finish();
        let original = SymbolicTensor::infer(expression).unwrap();
        let zero = original.with_checked_expression(Atom::Zero).unwrap();
        for source in [original, zero] {
            let mut calls = HashMap::<LibraryRep, usize>::new();
            let parametric = source
                .map_representation_dimensions(|rep, dimension| {
                    *calls.entry(rep).or_default() += 1;
                    if rep == (ColorAdjoint {}).into() {
                        Atom::var(CS.nc).pow(2) - 1
                    } else if rep == (ColorFundamental {}).into() {
                        Atom::var(CS.nc)
                    } else {
                        dimension.to_owned()
                    }
                })
                .unwrap();
            assert_eq!(calls[&(ColorAdjoint {}).into()], 1);
            assert_eq!(calls[&(ColorFundamental {}).into()], 1);
            assert_eq!(parametric.rank(), source.rank());
            assert_eq!(
                parametric.structure.open_positions(),
                source.structure.open_positions()
            );
            if !source.is_zero() {
                // Compact representation metadata must follow indexed and AUTO
                // ports alike; the old Atom color helper only matched indexed slots.
                // Dimension mapping encodes each representation before rebuilding
                // its Scalar container, which representation cooking leaves opaque.
                let expected_representation =
                    ColorAdjoint {}.to_symbolic([Atom::var(CS.nc).pow(2) - 1]);
                let expected_representation = CookSettings::indices()
                    .with_mode(CookMode::ReversibleEncoding)
                    .with_representation_payloads(true, true)
                    .try_cook(expected_representation.as_view())
                    .unwrap();
                let expected = FunctionBuilder::new(metadata)
                    .add_arg(expected_representation)
                    .finish();
                assert!(parametric.expression.contains(&expected));
            }
            assert!(
                parametric
                    .structure
                    .logical_slots()
                    .iter()
                    .filter(|slot| slot.rep().rep == (ColorAdjoint {}).into())
                    .all(|slot| {
                        crate::CookSettings::dimension_payload(
                            slot.rep().dim.to_symbolic().get_symbol().unwrap(),
                        ) == Some(Atom::var(CS.nc).pow(2) - 1)
                    })
            );
            let concrete = parametric
                .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))
                .unwrap();
            assert!(concrete.same_value(&source));
            assert!(
                source
                    .map_representation_dimensions(|_, _| {
                        FunctionBuilder::new(symbolica::symbol!("invalid_dimension_function"))
                            .add_arg(1)
                            .finish()
                    })
                    .is_err()
            );
        }
    }

    #[test]
    fn checked_symbol_renaming_preserves_swaps_zero_and_rejects_collisions() {
        crate::test_support::test_initialize();
        let a = symbolica::symbol!("checked_rename_a");
        let b = symbolica::symbol!("checked_rename_b");
        let dim = symbolica::symbol!("checked_rename_dimension");
        let other_dim = symbolica::symbol!("checked_rename_other_dimension");
        let rep = ExtendibleReps::MINKOWSKI.new_rep(dim);
        let source = SymbolicTensor::infer(
            FunctionBuilder::new(spenso::tensor_symbol!("checked_renaming_tensor"))
                .add_arg(rep.to_symbolic([Atom::var(a)]))
                .add_arg(rep.to_symbolic([Atom::var(b)]))
                .add_arg(rep.to_symbolic([]))
                .finish(),
        )
        .unwrap();
        for value in [
            source.clone(),
            source.with_checked_expression(Atom::Zero).unwrap(),
        ] {
            let renamed = value
                .rename_symbols(&HashMap::from([(a, b), (b, a), (dim, other_dim)]))
                .unwrap();
            let restored = renamed
                .rename_symbols(&HashMap::from([(a, b), (b, a), (other_dim, dim)]))
                .unwrap();
            assert!(restored.same_value(&value));
            assert_eq!(renamed.rank(), value.rank());
            assert_eq!(
                renamed.structure.open_positions(),
                value.structure.open_positions()
            );
            assert!(
                renamed
                    .structure
                    .logical_slots()
                    .iter()
                    .all(|slot| slot.rep().dim == Dimension::from(other_dim))
            );
            assert!(value.rename_symbols(&HashMap::from([(a, b)])).is_err());
            assert!(
                value
                    .rename_symbols(&HashMap::from([(a, b), (b, b)]))
                    .is_err()
            );
        }
    }

    #[test]
    fn symbol_renaming_preserves_compliant_callback_interface() {
        test_support::test_initialize();
        let source_index = symbol!("rename_callback_source");
        let target_index = symbol!("rename_callback_target");
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = tensor_symbol!(
            "rename_callback_tensor",
            norm = move |node, out| {
                observed.fetch_add(1, Ordering::Relaxed);
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(slot)) = function.iter().next()
                    && slot.iter().nth(1) == Some(Atom::var(target_index).as_view())
                {
                    **out = Atom::num(2)
                        * FunctionBuilder::new(tensor_symbol!("rename_compliant_result"))
                            .add_args(function.iter())
                            .finish();
                }
            }
        );
        let source = SymbolicTensor::infer(
            FunctionBuilder::new(head)
                .add_arg(
                    ExtendibleReps::MINKOWSKI
                        .new_rep(4)
                        .to_symbolic([Atom::var(source_index)]),
                )
                .finish(),
        )
        .unwrap();
        calls.store(0, Ordering::Relaxed);
        let renamed = source
            .rename_symbols(&HashMap::from([(source_index, target_index)]))
            .unwrap();
        assert_eq!(
            SymbolicTensor::infer(renamed.expression.clone())
                .unwrap()
                .structure,
            renamed.structure
        );
        assert_eq!(renamed.rank(), source.rank());
        assert_eq!(
            renamed.structure.logical_slots()[0].aind,
            PartialIndex::Explicit(AbstractIndex::Symbol(target_index.into()))
        );
        assert_eq!(calls.load(Ordering::Relaxed), 1);
        assert!(renamed.with_checked_expression(Atom::num(2)).is_err());
    }

    #[test]
    fn lorentz_dimension_merges_new_index_pairs_and_rejects_excess_occurrences() {
        let index = PartialIndex::Explicit(AbstractIndex::Normal(73311));
        let minkowski4 = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let minkowski6 = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(6));
        let remaining = ExtendibleReps::EUCLIDEAN
            .new_rep(Dimension::Concrete(3))
            .slot(PartialIndex::Explicit(AbstractIndex::Normal(73313)));
        let ports = [minkowski4.slot(index), remaining, minkowski6.slot(index)];
        let head = spenso::tensor_symbol!("dimension_collision_tensor");
        let expression = FunctionBuilder::new(head)
            .add_args(ports.map(composition::port_atom))
            .finish();
        let value = SymbolicTensor::<PartialStructure>::from_produced_parts(
            expression,
            PartialStructure::from_logical_slots(ports.into_iter().rev()),
        )
        .unwrap();
        assert_eq!(value.rank(), 3);
        let result = value
            .with_lorentz_dimension(Dimension::Concrete(6))
            .unwrap();
        assert_eq!(result.structure.logical_slots(), vec![remaining]);
        let ports = [
            minkowski4.slot(index),
            minkowski6.slot(index),
            minkowski6.slot(index),
        ];
        let expression = FunctionBuilder::new(head)
            .add_args(ports.map(composition::port_atom))
            .finish();
        let value = SymbolicTensor::<PartialStructure>::from_produced_parts(
            expression,
            PartialStructure::from_logical_slots(ports),
        )
        .unwrap()
        .indexed(&HashMap::new())
        .unwrap();
        assert_eq!(value.rank(), 1);
        assert!(
            value
                .with_lorentz_dimension(Dimension::Concrete(6))
                .is_err()
        );
    }

    #[test]
    fn dimension_mapping_preserves_compliant_callback_interface() {
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let head = tensor_symbol!(
            "dimension_result_normalizer",
            norm = move |node, output| {
                observed.lock().unwrap().push(node.to_owned());
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(slot)) = function.iter().last()
                    && slot.iter().next() == Some(Atom::num(6).as_view())
                {
                    let mode = function.iter().next().unwrap();
                    if mode == Atom::Zero.as_view() || mode == Atom::one().as_view() {
                        **output = if mode.is_zero() {
                            Atom::Zero
                        } else {
                            Atom::num(2)
                                * FunctionBuilder::new(tensor_symbol!("dimension_compliant_result"))
                                    .add_args(function.iter())
                                    .finish()
                        };
                    }
                }
            }
        );
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        for index in [
            PartialIndex::open(0),
            PartialIndex::Explicit(AbstractIndex::Normal(73321)),
        ] {
            for mode in [0, 1, 2] {
                let port = representation.slot(index);
                let expression = FunctionBuilder::new(head)
                    .add_arg(mode)
                    .add_arg(composition::port_atom(port))
                    .finish();
                let value = SymbolicTensor::<PartialStructure>::new(
                    expression,
                    PartialStructure::from_logical_slots([port]),
                )
                .unwrap();
                calls.lock().unwrap().clear();
                let unchanged = value
                    .expression
                    .with_lorentz_dimension(Atom::num(4).as_view());
                let unchanged_calls = mem::take(&mut *calls.lock().unwrap());
                let mut mapped_pairs = Vec::new();
                let actual = value
                    .map_representation_dimensions(|representation, dimension| {
                        mapped_pairs.push((representation, dimension.to_owned()));
                        dimension.to_owned()
                    })
                    .unwrap();
                assert_eq!(actual.expression, unchanged);
                assert_eq!(actual.structure, value.structure);
                assert_eq!(*calls.lock().unwrap(), unchanged_calls);
                assert_eq!(mapped_pairs, vec![(representation.rep, Atom::num(4))]);

                calls.lock().unwrap().clear();
                let expected = value
                    .expression
                    .with_lorentz_dimension(Atom::num(6).as_view());
                let expected_calls = mem::take(&mut *calls.lock().unwrap());
                assert!(!expected_calls.is_empty());
                let result = value.with_lorentz_dimension(Dimension::Concrete(6));
                assert_eq!(*calls.lock().unwrap(), expected_calls);
                let result = result.unwrap();
                assert_eq!(result.expression, expected);
                if !result.is_zero() {
                    assert_eq!(
                        SymbolicTensor::infer(result.expression.clone())
                            .unwrap()
                            .structure,
                        result.structure
                    );
                }
                assert_eq!(
                    result.structure.logical_slots(),
                    vec![
                        ExtendibleReps::MINKOWSKI
                            .new_rep(Dimension::Concrete(6))
                            .slot(index)
                    ]
                );
            }
        }

        let promoted_tensor = Atom::num(7);
        let observed = Arc::clone(&calls);
        let scalar = symbol!(
            "dimension_scalar_callback"; Scalar;
            norm = move |node, output| {
                observed.lock().unwrap().push(node.to_owned());
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(slot)) = function.iter().next()
                    && slot.iter().next() == Some(Atom::num(6).as_view())
                {
                    **output = promoted_tensor.clone();
                }
            }
        );
        let expression = FunctionBuilder::new(scalar)
            .add_arg(representation.to_symbolic([]))
            .finish();
        assert!(!InterfaceInference::has_structured_syntax(
            expression.as_view()
        ));
        let value = SymbolicTensor::<PartialStructure>::new(
            expression,
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        calls.lock().unwrap().clear();
        let expected = value
            .expression
            .with_lorentz_dimension(Atom::num(6).as_view());
        assert_eq!(expected, Atom::num(7));
        let expected_calls = mem::take(&mut *calls.lock().unwrap());
        let mapped = value
            .with_lorentz_dimension(Dimension::Concrete(6))
            .unwrap();
        assert_eq!(mapped.expression, Atom::num(7));
        assert!(mapped.is_scalar());
        assert_eq!(*calls.lock().unwrap(), expected_calls);
    }

    #[test]
    fn cooking_preserves_logical_order_and_rejects_invalid_source() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let named = symbol!("shared_cooking_index", tags = [&SPENSO_TAG.index]);
        let slots = [73101, 73103].map(|owner| {
            representation.slot(PartialIndex::Explicit(AbstractIndex::Named(
                named.into(),
                owner,
                2,
            )))
        });
        let expression = FunctionBuilder::new(tensor_symbol!("shared_cooking_tensor"))
            .add_args(slots.map(composition::port_atom))
            .finish();
        let original = SymbolicTensor::<PartialStructure>::from_declared_parts(
            expression.clone(),
            PartialStructure::from_logical_slots(slots.into_iter().rev()),
        )
        .unwrap();
        let settings = CookSettings::indices();
        let expected = slots
            .into_iter()
            .rev()
            .map(composition::port_atom)
            .collect::<Vec<_>>();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let cooked = original.with_cooked_indices(&settings).unwrap();
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        assert_eq!(cooked.expression, expression);
        assert_eq!(
            cooked
                .structure
                .logical_slots()
                .into_iter()
                .map(composition::port_atom)
                .collect::<Vec<_>>(),
            expected
        );

        let unknown = FunctionBuilder::new(symbol!("cooking_callback_unknown_index"))
            .add_arg(2)
            .finish();
        let raw_slot = representation.to_symbolic([unknown]);
        let target = settings.try_cook_indices(raw_slot.as_view()).unwrap();
        assert_ne!(target, raw_slot);
        let slot = Slot::<LibraryRep, AbstractIndex>::try_from(target.as_view()).unwrap();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let callback = tensor_symbol!(
            "cooking_callback_strips_ports",
            norm = move |node, output| {
                observed.fetch_add(1, Ordering::Relaxed);
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let expression = FunctionBuilder::new(callback).add_arg(raw_slot).finish();
        let structure = PartialStructure::from_logical_slots([slot
            .rep()
            .slot(PartialIndex::Explicit(slot.aind()))]);
        calls.store(0, Ordering::Relaxed);
        let raw_error = SymbolicTensor::new(expression.clone(), structure.clone()).unwrap_err();
        assert!(matches!(raw_error, TensorInferenceError::Invalid(message)
            if message.contains("cook nested indices")));
        assert_eq!(calls.load(Ordering::Relaxed), 0);
        let cooked = settings.try_cook_indices(expression.as_view()).unwrap();
        assert!(cooked.is_one());
        assert_eq!(calls.load(Ordering::Relaxed), 1);
        let cooked_error = SymbolicTensor::from_declared_parts(cooked, structure).unwrap_err();
        assert!(
            matches!(cooked_error, TensorInferenceError::Invalid(message)
            if message.contains("compatible tensor interface"))
        );
        assert_eq!(calls.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn typed_scoping_preserves_compliant_callback_interface_and_ranked_zero() {
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let scope = symbol!("shared_scope_test");
        let index = AbstractIndex::Normal(73201);
        let slot = rep.slot(PartialIndex::Explicit(index));
        let target = rep.slot::<AbstractIndex, _>(index.scoped(scope)).to_atom();
        let callback = tensor_symbol!(
            "scope_callback_strips_ports",
            norm = move |node, output| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|arg| arg == target.as_view())
                {
                    **output = Atom::num(2)
                        * FunctionBuilder::new(tensor_symbol!("scope_compliant_result"))
                            .add_args(function.iter())
                            .finish();
                }
            }
        );
        let value = SymbolicTensor::<PartialStructure>::new(
            FunctionBuilder::new(callback)
                .add_arg(composition::port_atom(slot))
                .finish(),
            PartialStructure::from_logical_slots([slot]),
        )
        .unwrap();
        let scoped = value.wrap_indices(scope, false).unwrap();
        assert_eq!(
            SymbolicTensor::infer(scoped.expression.clone())
                .unwrap()
                .structure,
            scoped.structure
        );
        assert_eq!(
            scoped.structure.logical_slots(),
            [rep.slot(PartialIndex::Explicit(index.scoped(scope)))]
        );
        assert_eq!(
            value.wrap_indices(scope, true).unwrap().expression(),
            value.expression()
        );
        let zero = SymbolicTensor::<PartialStructure>::new(
            Atom::Zero,
            PartialStructure::from_logical_slots([slot, rep.slot(PartialIndex::open(1))]),
        )
        .unwrap();
        let scoped = zero.wrap_indices(scope, false).unwrap();
        assert!(scoped.expression().is_zero());
        assert_eq!(
            scoped.structure.logical_slots()[0].aind,
            PartialIndex::Explicit(index.scoped(scope))
        );
        assert_eq!(scoped.structure.open_positions(), vec![1]);
        assert_eq!(
            zero.wrap_indices(scope, true).unwrap().structure,
            zero.structure
        );
    }

    #[test]
    fn checked_addition_owns_canonical_axes_and_typed_zero_without_reinference() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let slots = [73111, 73113]
            .map(|index| representation.slot(PartialIndex::Explicit(AbstractIndex::Normal(index))));
        let expression = FunctionBuilder::new(tensor_symbol!("shared_additive_tensor"))
            .add_args(slots.map(composition::port_atom))
            .finish();
        let left = SymbolicTensor::<PartialStructure>::from_declared_parts(
            expression.clone(),
            PartialStructure::from_logical_slots(slots.into_iter().rev()),
        )
        .unwrap();
        let right = SymbolicTensor::<PartialStructure>::from_declared_parts(
            expression,
            PartialStructure::from_logical_slots(slots),
        )
        .unwrap();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let sum = left.try_add(&right).unwrap();
        let zero = left.try_sub(&right).unwrap();
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        assert_eq!(
            sum.structure.logical_slots(),
            left.structure.canonical().external_structure()
        );
        assert_eq!(
            zero.structure.logical_slots(),
            left.structure.canonical().external_structure()
        );
        assert!(zero.expression.is_zero());
        let scalar = SymbolicTensor::<PartialStructure>::new(
            Atom::one(),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        assert!(left.try_add(&scalar).is_err());
        assert!(left.try_sub(&scalar).is_err());
    }

    #[test]
    fn metric_result_validation_rejects_callback_rank_loss_and_retains_typed_zero() {
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let a = representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(72001))
            .to_atom();
        let b = representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(72003))
            .to_atom();
        let target = b.clone();
        let tensor = spenso::tensor_symbol!(
            "metric_inference_callback_rank_loss",
            norm = move |node, output| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let atom = FunctionBuilder::new(ETS.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish()
            * FunctionBuilder::new(tensor).add_arg(a).finish();
        let value = SymbolicTensor::<PartialStructure>::infer(atom.clone()).unwrap();
        assert_eq!(value.rank(), 1);
        SymbolicTensor::<PartialStructure>::validate_interface(&atom, &value.structure).unwrap();

        let contracted =
            atom.schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors());
        assert_eq!(contracted, Atom::one());
        assert!(
            SymbolicTensor::<PartialStructure>::validate_interface(&contracted, &value.structure)
                .is_err()
        );
        assert!(value.with_checked_expression(contracted).is_err());
        SymbolicTensor::<PartialStructure>::validate_interface(&Atom::Zero, &value.structure)
            .unwrap();
    }

    #[test]
    fn metric_gamma_rewrites_preserve_auto_ports_and_logical_order() {
        // This test selects family prerequisites only, preserving its explicit output scope.

        use crate::shorthands::schoonschip::Schoonschip;

        crate::test_support::test_initialize();
        let mink = LibraryRep::from(Minkowski {}).new_rep(4);
        let [mu, nu] = [74801, 74803].map(AbstractIndex::Normal);
        let metric = SymbolicTensor::<PartialStructure>::infer(
            FunctionBuilder::new(ETS.metric)
                .add_args([mu, nu].map(|index| mink.slot::<AbstractIndex, _>(index).to_atom()))
                .finish(),
        )
        .unwrap();
        for explicit_spin in [false, true] {
            let mut indices = HashMap::from([(2, nu)]);
            if explicit_spin {
                indices.extend([
                    (0, AbstractIndex::Normal(74805)),
                    (1, AbstractIndex::Normal(74807)),
                ]);
            }
            let gamma = SymbolicTensor::<PartialStructure>::from_signature(
                &AGS.gamma_strct::<AbstractIndex>(4),
            )
            .unwrap()
            .reindex_interface_ports(&indices)
            .unwrap();
            let expected = gamma
                .reindex_interface_ports(&HashMap::from([(2, mu)]))
                .unwrap();
            for source in [
                metric.multiply(&gamma).unwrap(),
                gamma.multiply(&metric).unwrap(),
            ] {
                let original = source.expression.clone();
                for (rewritten, expected) in [
                    (
                        source.expression.schoonschip_with_settings(
                            &SchoonschipSettings::default().without_rank1_tensors(),
                        ),
                        expected.expression.clone(),
                    ),
                    (source.expression.schoonschip(), expected.expression.clone()),
                    (
                        crate::tensor::SymbolicTensor::infer(
                            (source.expression).as_atom_view().to_owned(),
                        )
                        .unwrap()
                        .simplify_algebra(&crate::tensor::AlgebraSettings {
                            contract: crate::tensor::AlgebraContraction::None,
                            gamma: Some(crate::dirac::GammaSimplifySettings::default()),
                            epsilon: true,
                            ..Default::default()
                        })
                        .unwrap()
                        .contract(crate::tensor::ContractSettings {
                            collect_chains: false,
                            collect_traces: false,
                            ..Default::default()
                        })
                        .unwrap()
                        .into_expression(),
                        crate::tensor::SymbolicTensor::infer(
                            (expected.expression).as_atom_view().to_owned(),
                        )
                        .unwrap()
                        .simplify_algebra(&crate::tensor::AlgebraSettings {
                            contract: crate::tensor::AlgebraContraction::None,
                            gamma: Some(crate::dirac::GammaSimplifySettings::default()),
                            epsilon: true,
                            ..Default::default()
                        })
                        .unwrap()
                        .contract(crate::tensor::ContractSettings {
                            collect_chains: false,
                            collect_traces: false,
                            ..Default::default()
                        })
                        .unwrap()
                        .into_expression(),
                    ),
                ] {
                    assert_eq!(rewritten, expected);
                    let result = source.with_rewritten_expression(rewritten.clone()).unwrap();
                    assert_eq!(result.expression, rewritten);
                    assert_eq!(
                        result.structure.logical_slots(),
                        source.structure.logical_slots()
                    );
                    assert_eq!(
                        result.structure.open_positions().len(),
                        if explicit_spin { 0 } else { 2 }
                    );
                    assert_eq!(result.with_rewritten_expression(rewritten).unwrap(), result);
                }
                assert_eq!(source.expression, original);
                let zero = source.with_rewritten_expression(Atom::Zero).unwrap();
                assert_eq!(
                    zero.structure.logical_slots(),
                    source.structure.logical_slots()
                );
            }
        }
    }

    #[test]
    fn rewritten_auto_ports_require_the_original_encoded_owner_and_axis() {
        crate::test_support::test_initialize();
        let mink = LibraryRep::from(Minkowski {}).new_rep(4);
        let spin = LibraryRep::from(Bispinor {}).new_rep(4);
        let owner = AbstractIndex::fresh_open_owner();
        let ends = [0, 1].map(|axis| AbstractIndex::Open { owner, axis });
        let mu = AbstractIndex::Normal(74811);
        let nu = AbstractIndex::Normal(74813);
        let gamma = |ends: [AbstractIndex; 2], index| {
            FunctionBuilder::new(AGS.gamma)
                .add_args(ends.map(|index| spin.slot::<AbstractIndex, _>(index).to_atom()))
                .add_arg(mink.slot::<AbstractIndex, _>(index).to_atom())
                .finish()
        };
        let metric = FunctionBuilder::new(ETS.metric)
            .add_args([mu, nu].map(|index| mink.slot::<AbstractIndex, _>(index).to_atom()))
            .finish();
        let logical = PartialStructure::from_logical_slots([
            mink.slot(PartialIndex::Explicit(mu)),
            spin.slot(PartialIndex::open(0)),
            spin.slot(PartialIndex::open(1)),
        ]);
        let source =
            SymbolicTensor::from_produced_parts(metric * gamma(ends, nu), logical).unwrap();
        let result = gamma(ends, mu);
        // Both boundaries can retain the public order when encoded owners prove
        // that moving the explicit Lorentz port past AUTO axes is harmless.
        for finish in [
            SymbolicTensor::<PartialStructure>::with_checked_expression,
            SymbolicTensor::<PartialStructure>::with_rewritten_expression,
        ] {
            let accepted = finish(&source, result.clone()).unwrap();
            assert_eq!(accepted.expression, result);
            assert_eq!(
                accepted.structure.logical_slots(),
                source.structure.logical_slots()
            );
        }
        for changed in [
            AbstractIndex::Open {
                owner: owner + 1,
                axis: 0,
            },
            AbstractIndex::Open { owner, axis: 2 },
        ] {
            assert!(
                source
                    .with_checked_expression(gamma([changed, ends[1]], mu))
                    .is_err()
            );
        }
        let wrong_branch = gamma(ends, AbstractIndex::Normal(74815));
        assert!(
            source
                .with_checked_expression(gamma(ends, mu) + wrong_branch)
                .is_err()
        );
        let reordered = PartialStructure::from_logical_slots([
            spin.slot(PartialIndex::open(0)),
            spin.slot(PartialIndex::open(1)),
            mink.slot(PartialIndex::Explicit(mu)),
        ]);
        assert!(!InterfaceInference::additive_interfaces_match(
            &source.structure,
            &reordered
        ));
        let anonymous = FunctionBuilder::new(spenso::tensor_symbol!("rewrite_anonymous_ports"))
            .add_arg(spin.to_symbolic([]))
            .add_arg(mink.slot::<AbstractIndex, _>(mu).to_atom())
            .finish();
        let invalid = SymbolicTensor::new(
            anonymous,
            PartialStructure::from_logical_slots([
                mink.slot(PartialIndex::Explicit(mu)),
                spin.slot(PartialIndex::open(0)),
            ]),
        )
        .unwrap_err();
        assert!(matches!(invalid, TensorInferenceError::Invalid(message)
            if message.contains("compatible tensor interface")));
    }

    #[test]
    fn encoded_rewrite_observation_rejects_callback_rank_loss_without_replay() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };

        crate::test_support::test_initialize();
        let mink = LibraryRep::from(Minkowski {}).new_rep(4);
        let spin = LibraryRep::from(Bispinor {}).new_rep(4);
        let [mu, nu] = [74821, 74823].map(|i| {
            mink.slot::<AbstractIndex, _>(AbstractIndex::Normal(i))
                .to_atom()
        });
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let target = mu.clone();
        let head = spenso::tensor_symbol!(
            "auto_metric_callback_rank_loss",
            norm = move |node, out| {
                observed.fetch_add(1, Ordering::Relaxed);
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **out = Atom::one();
                }
            }
        );
        let owner = AbstractIndex::fresh_open_owner();
        let tensor = SymbolicTensor::<PartialStructure>::infer(
            FunctionBuilder::new(head)
                .add_args([0, 1].map(|axis| {
                    spin.slot::<AbstractIndex, _>(AbstractIndex::Open { owner, axis })
                        .to_atom()
                }))
                .add_arg(&nu)
                .finish(),
        )
        .unwrap();
        let metric = SymbolicTensor::<PartialStructure>::infer(
            FunctionBuilder::new(ETS.metric)
                .add_arg(&mu)
                .add_arg(&nu)
                .finish(),
        )
        .unwrap();
        let source = metric.multiply(&tensor).unwrap();
        let result = source
            .expression
            .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors());
        assert_eq!(result, Atom::one());
        let before = calls.load(Ordering::Relaxed);
        assert!(before > 0);
        assert!(source.with_checked_expression(result).is_err());
        assert_eq!(calls.load(Ordering::Relaxed), before);
    }

    #[test]
    fn terminal_trace_rewrites_retain_admitted_structure_and_logical_order() {
        crate::test_support::test_initialize();
        let spin = Bispinor {}.new_rep(4).to_symbolic([]);
        let spectator = symbolica::parse_lit!((trace_interface_x + trace_interface_y) ^ 3);
        for dimension in [
            Dimension::Concrete(4),
            Dimension::from(symbolica::symbol!("trace_interface_D")),
        ] {
            let minkowski = Minkowski {}.new_rep(dimension);
            let indices = (0..8)
                .map(|position| {
                    minkowski
                        .slot::<AbstractIndex, _>(if position < 4 {
                            AbstractIndex::Normal(73501 + position)
                        } else {
                            AbstractIndex::Dummy(73501 + position)
                        })
                        .to_atom()
                })
                .collect::<Vec<_>>();
            let vector = spenso::vector_symbol!("trace_interface_vector");
            let compact = FunctionBuilder::new(vector)
                .add_arg(minkowski.to_symbolic([]))
                .finish();
            let distinct_vectors = (0..8)
                .map(|position| {
                    FunctionBuilder::new(vector)
                        .add_arg(position)
                        .add_arg(minkowski.to_symbolic([]))
                        .finish()
                })
                .collect::<Vec<_>>();
            let words = [
                indices.clone(),
                distinct_vectors,
                vec![
                    indices[0].clone(),
                    compact.clone(),
                    indices[0].clone(),
                    compact.clone(),
                ],
                vec![compact.clone(), compact.clone(), compact.clone(), compact],
            ];
            for word in words {
                let ordinary =
                    shadowing::trace(&spin, word.iter().map(|index| crate::gamma!(index)));
                let mut expressions = vec![ordinary.clone(), &spectator * ordinary];
                if dimension == Dimension::Concrete(4) && word == indices {
                    let axial = shadowing::trace(
                        &spin,
                        std::iter::once(crate::gamma5!())
                            .chain(word.iter().map(|index| crate::gamma!(index))),
                    );
                    expressions.extend([axial.clone(), &spectator * axial]);
                }
                for expression in expressions {
                    let mut value = SymbolicTensor::<PartialStructure>::infer(expression).unwrap();
                    value.structure = PartialStructure::from_logical_slots(
                        value.structure.logical_slots().into_iter().rev(),
                    );
                    let expected = crate::tensor::SymbolicTensor::infer(
                        (value.expression).as_atom_view().to_owned(),
                    )
                    .unwrap()
                    .simplify_algebra(&crate::tensor::AlgebraSettings {
                        gamma: Some(crate::dirac::GammaSimplifySettings::default()),
                        epsilon: true,
                        ..Default::default()
                    })
                    .unwrap()
                    .contract(crate::tensor::ContractSettings {
                        collect_chains: false,
                        collect_traces: false,
                        ..Default::default()
                    })
                    .unwrap()
                    .into_expression();
                    assert_ne!(expected, value.expression);
                    if word.len() == 8 {
                        assert!(
                            expected.as_view().get_byte_size()
                                > value.expression.as_view().get_byte_size().saturating_mul(2)
                        );
                    }
                    INFERENCE_CALLS.with(|count| count.set(0));
                    INFERENCE_VISITS.with(|count| count.set(0));
                    let result = value.with_rewritten_expression(expected.clone()).unwrap();
                    assert_eq!(result.expression, expected);
                    assert_eq!(result.structure, value.structure);
                    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
                    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
                    SymbolicTensor::<PartialStructure>::validate_interface(
                        &result.expression,
                        &result.structure,
                    )
                    .unwrap();
                }
            }
        }
    }

    #[test]
    fn terminal_trace_rejects_user_callback_rank_loss() {
        use std::sync::{Arc, Mutex};

        crate::test_support::test_initialize();
        let scalar = Atom::add_many((0..64).map(|index| {
            Atom::var(symbolica::symbol!(format!(
                "trace_interface_callback_scalar_{index}"
            )))
        }));
        let callback_scalar = scalar.clone();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let recorded = Arc::clone(&calls);
        let vector = spenso::vector_symbol!(
            "trace_interface_callback_vector",
            norm = move |node, out| {
                recorded.lock().unwrap().push(node.to_owned());
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(slot)) = function.iter().last()
                    && slot.get_nargs() == 2
                {
                    **out = callback_scalar.clone();
                }
            }
        );
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(73531);
        let explicit = minkowski.slot::<AbstractIndex, _>(index).to_atom();
        let compact = FunctionBuilder::new(vector)
            .add_arg(minkowski.to_symbolic([]))
            .finish();
        let source = shadowing::trace(
            Bispinor {}.new_rep(4).to_symbolic([]),
            [&explicit, &compact].map(|index| crate::gamma!(index)),
        );
        let value =
            SymbolicTensor::<PartialStructure>::new(source, explicit_interface(minkowski, index))
                .unwrap();
        calls.lock().unwrap().clear();
        let settings = crate::dirac::GammaSimplifySettings::default();
        let rewritten = crate::dirac::DiracSimplifier::new(&settings)
            .evaluate_terminal_trace::<false>(value.expression.as_view())
            .map(|(expression, _)| expression)
            .unwrap();
        assert_eq!(rewritten, Atom::num(4) * scalar);
        assert!(
            rewritten.as_view().get_byte_size()
                > value.expression.as_view().get_byte_size().saturating_mul(2)
        );
        let before_validation = calls.lock().unwrap().clone();
        assert!(!before_validation.is_empty());
        assert!(value.with_checked_expression(rewritten).is_err());
        assert_eq!(*calls.lock().unwrap(), before_validation);
        assert!(
            value
                .simplify_algebra(&crate::tensor::AlgebraSettings {
                    gamma: Some(settings),
                    epsilon: settings.output == crate::dirac::GammaOutput::Reduced,
                    ..Default::default()
                })
                .is_err()
        );
    }

    #[test]
    fn proven_algebra_metric_and_zero_results_retain_interfaces_without_reinference() {
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let [a, b] = [72011, 72013].map(|index| {
            representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let tensor = tensor_symbol!("shared_reuse_tensor");
        let vector = |slot: &Atom| FunctionBuilder::new(tensor).add_arg(slot).finish();
        let metric = FunctionBuilder::new(ETS.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish();
        let original = SymbolicTensor::<PartialStructure>::infer(metric * vector(&a)).unwrap();
        let rewritten = original
            .expression
            .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors());
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let result = original
            .with_rewritten_expression(rewritten.clone())
            .unwrap();
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        assert_eq!(result.expression, vector(&b));
        assert_eq!(
            result.structure.logical_slots(),
            original.structure.logical_slots()
        );

        let x = Atom::var(symbol!("shared_reuse_x"));
        let y = Atom::var(symbol!("shared_reuse_y"));
        let algebra = SymbolicTensor::<PartialStructure>::infer((x + y) * vector(&a)).unwrap();
        let expanded = algebra.expression.expand();
        assert_ne!(expanded, algebra.expression);
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let expanded = algebra.with_algebra_result(expanded).unwrap();
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        assert_eq!(
            expanded.structure.logical_slots(),
            algebra.structure.logical_slots()
        );

        let unresolved = SymbolicTensor::<PartialStructure>::from_produced_parts(
            FunctionBuilder::new(tensor)
                .add_arg(representation.to_symbolic([]))
                .finish(),
            PartialStructure::from_logical_slots([representation.slot(PartialIndex::open(0))]),
        )
        .unwrap();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        for result in [
            unresolved.with_algebra_result(Atom::Zero).unwrap(),
            unresolved.with_rewritten_expression(Atom::Zero).unwrap(),
            unresolved
                .with_rewritten_expression(unresolved.expression.clone())
                .unwrap(),
        ] {
            assert_eq!(
                result.structure.logical_slots(),
                unresolved.structure.logical_slots()
            );
        }
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn unchanged_results_preserve_layout_dispatch_flags_and_callback_state() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };

        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = spenso::tensor_symbol!(
            "unchanged_result_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let ports = [
            ExtendibleReps::MINKOWSKI
                .new_rep(Dimension::Concrete(4))
                .slot(PartialIndex::Explicit(AbstractIndex::Normal(72017))),
            ExtendibleReps::EUCLIDEAN
                .new_rep(Dimension::Concrete(3))
                .slot(PartialIndex::open(0)),
        ];
        let expression = FunctionBuilder::new(head)
            .add_arg(17)
            .add_args(ports.map(composition::port_atom))
            .finish();
        for expression in [expression, Atom::Zero] {
            for (is_metric, is_composite) in
                [(false, false), (false, true), (true, false), (true, true)]
            {
                let mut value = SymbolicTensor::from_declared_parts(
                    expression.clone(),
                    PartialStructure::from_logical_slots(ports),
                )
                .unwrap();
                value.is_metric = is_metric;
                value.is_composite = is_composite;
                for finish in [
                    SymbolicTensor::<PartialStructure>::with_checked_expression,
                    SymbolicTensor::<PartialStructure>::with_algebra_result,
                    SymbolicTensor::<PartialStructure>::with_rewritten_expression,
                    |source: &SymbolicTensor<PartialStructure>, result: Atom| {
                        source
                            .with_transformed_expression(result)
                            .map(|(value, contracted)| {
                                assert!(!contracted);
                                value
                            })
                    },
                ] {
                    calls.store(0, Ordering::Relaxed);
                    INFERENCE_CALLS.with(|count| count.set(0));
                    INFERENCE_VISITS.with(|count| count.set(0));
                    let result = finish(&value, value.expression.clone()).unwrap();
                    assert_eq!(result, value);
                    assert_eq!(calls.load(Ordering::Relaxed), 0);
                    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
                    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
                }
            }
        }
    }

    #[test]
    fn arbitrary_transform_reports_raw_contractions_without_losing_typed_zero() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbolica::symbol!("transform_policy_x"));
        let scalar = SymbolicTensor::from_produced_parts(
            &x + Atom::one(),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        let slot = ExtendibleReps::EUCLIDEAN
            .new_rep(Dimension::Concrete(3))
            .slot(PartialIndex::Explicit(AbstractIndex::Normal(73201)));
        let head = spenso::tensor_symbol!("transform_policy_tensor");
        let diagonal = FunctionBuilder::new(head)
            .add_arg(composition::port_atom(slot))
            .add_arg(composition::port_atom(slot))
            .finish();
        let (value, contracted) = scalar.with_transformed_expression(diagonal).unwrap();
        assert_eq!(value.rank(), scalar.rank());
        assert!(
            contracted,
            "a data descriptor must be discarded even at equal final rank"
        );
        let (value, contracted) = scalar
            .with_transformed_expression(&x + Atom::num(2))
            .unwrap();
        assert_eq!(value.rank(), 0);
        assert!(!contracted);

        let tensor = SymbolicTensor::from_produced_parts(
            FunctionBuilder::new(head)
                .add_arg(composition::port_atom(slot))
                .finish(),
            PartialStructure::from_logical_slots([slot]),
        )
        .unwrap();
        let (zero, contracted) = tensor.with_transformed_expression(Atom::Zero).unwrap();
        assert_eq!(zero.structure, tensor.structure);
        assert!(zero.expression.is_zero());
        assert!(!contracted);
        assert!(tensor.with_transformed_expression(Atom::one()).is_err());
    }

    #[test]
    fn reused_reusable_interfaces_preserve_interface_and_multiplicity_validation() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let [a, b, c] = [71, 73, 79].map(|index| {
            representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let metric = |left: &Atom, right: &Atom| {
            FunctionBuilder::new(ETS.metric)
                .add_arg(left)
                .add_arg(right)
                .finish()
        };
        let ab = metric(&a, &b);
        let ac = metric(&a, &c);
        let mut inference = InterfaceInference::default();
        let expected = inference.infer(&ab).unwrap();
        assert_eq!(
            inference.infer(&ab).unwrap().logical_slots(),
            expected.logical_slots()
        );
        assert!(
            inference
                .reusable_interfaces
                .contains_key(ab.as_view().get_data())
        );

        let contracted = inference.infer(&(ab.as_ref() * ac.as_ref())).unwrap();
        assert_eq!(contracted.canonical().order(), 2);
        for invalid in [
            ab.as_ref() + ac.as_ref(),
            ab.as_ref() + Atom::num(1),
            ab.as_ref().pow(Atom::num(3)),
        ] {
            assert!(inference.infer(&invalid).is_err());
        }

        // A leaf's own repeated slots must remain visible to enclosing products.
        let diagonal = FunctionBuilder::new(ETS.flat)
            .add_arg(&a)
            .add_arg(&a)
            .finish();
        assert_eq!(inference.infer(&diagonal).unwrap().canonical().order(), 2);
        let vector = FunctionBuilder::new(spenso::vector_symbol!("inference_reuse_vector"))
            .add_arg(&a)
            .finish();
        assert!(inference.infer(&(diagonal * vector)).is_err());

        let flat = FunctionBuilder::new(ETS.flat)
            .add_arg(&a)
            .add_arg(&b)
            .finish();
        let sum = ab.as_ref() * ac.as_ref() + flat * ac.as_ref();
        let actual = InterfaceInference::default()
            .infer_validated(sum.as_view())
            .unwrap();
        assert!(InterfaceInference::additive_interfaces_match(
            &contracted,
            &actual
        ));
    }

    #[test]
    fn interface_reuse_excludes_unresolved_builtin_ports() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let compact = representation.to_symbolic([]);
        let implicit = FunctionBuilder::new(ETS.metric)
            .add_arg(&compact)
            .add_arg(&compact)
            .finish();
        let materialized = FunctionBuilder::new(ETS.metric)
            .add_args((0..2).map(|axis| {
                representation
                    .slot::<AbstractIndex, _>(AbstractIndex::Open { owner: 81, axis })
                    .to_atom()
            }))
            .finish();
        let mut inference = InterfaceInference::default();
        for (atom, rank) in [(implicit, 2), (materialized, 2)] {
            for _ in 0..2 {
                let interface = inference.infer(&atom).unwrap();
                assert_eq!(interface.open_positions(), (0..rank).collect::<Vec<_>>());
                assert!(inference.reusable_interfaces.is_empty());
            }
        }
    }

    #[test]
    fn direct_leaf_inference_preserves_open_occurrences_and_logical_order() {
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let leaf = FunctionBuilder::new(spenso::tensor_symbol!("direct_open_leaf"))
            .add_arg(7)
            .add_arg(minkowski.to_symbolic([]))
            .add_arg(euclidean.to_symbolic([]))
            .finish();
        let mut inference = InterfaceInference::default();
        let interface = inference.infer_validated(leaf.as_view()).unwrap();
        assert_eq!(
            interface
                .logical_slots()
                .iter()
                .map(IsAbstractSlot::rep)
                .collect::<Vec<_>>(),
            vec![minkowski, euclidean],
        );
        assert!(
            inference
                .reusable_interfaces
                .contains_key(leaf.as_view().get_data())
        );
        let repeated = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(&leaf)
            .add_arg(&leaf)
            .finish();
        let interface = inference.infer_validated(repeated.as_view()).unwrap();
        assert_eq!(interface.open_positions(), vec![0, 1, 2, 3]);
        assert_eq!(
            interface
                .logical_slots()
                .into_iter()
                .map(|slot| slot.aind)
                .collect::<Vec<_>>(),
            (0..4).map(PartialIndex::open).collect::<Vec<_>>(),
        );
        let incompatible = FunctionBuilder::new(spenso::tensor_symbol!("direct_other_leaf"))
            .add_arg(euclidean.to_symbolic([]))
            .add_arg(minkowski.to_symbolic([]))
            .finish();
        assert!(
            inference
                .infer_validated((&leaf + incompatible).as_view())
                .is_err()
        );

        let p = FunctionBuilder::new(spenso::vector_symbol!("direct_compact_p"))
            .add_arg(minkowski.to_symbolic([]))
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("direct_compact_q"))
            .add_arg(minkowski.to_symbolic([]))
            .finish();
        assert!(
            inference
                .infer_validated(ETS.metric(p, q).as_view())
                .unwrap()
                .canonical()
                .is_scalar()
        );
    }

    #[test]
    fn direct_leaf_inference_retains_materialization_callback_semantics() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let callback = spenso::tensor_symbol!(
            "materialized_leaf_callback",
            norm = |node, output| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| {
                        Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok()
                    })
                {
                    **output = Atom::Zero;
                }
            }
        );
        let leaf = FunctionBuilder::new(callback)
            .add_arg(representation.to_symbolic([]))
            .finish();
        assert!(matches!(leaf.as_view(), AtomView::Fun(_)));
        let mut inference = InterfaceInference::default();
        assert!(
            inference
                .infer_validated(leaf.as_view())
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(inference.reusable_interfaces.is_empty());
        let value = SymbolicTensor::<PartialStructure>::infer(leaf).unwrap();
        // Intake observes the already constructed leaf. A later callback may
        // return zero, whose nominal port remains owned by this carrier.
        assert_eq!(value.rank(), 1);
        assert!(value.proofs.is_admitted());
        assert!(value.established_interface_is_valid());
    }

    #[test]
    fn repeated_scalar_dot_inference_preserves_external_ports_and_errors() {
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let p = FunctionBuilder::new(spenso::vector_symbol!("cached_dot_p"))
            .add_arg(representation.to_symbolic([]))
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("cached_dot_q"))
            .add_arg(representation.to_symbolic([]))
            .finish();
        let incompatible = FunctionBuilder::new(spenso::vector_symbol!("cached_dot_q"))
            .add_arg(ExtendibleReps::EUCLIDEAN.new_rep(4).to_symbolic([]))
            .finish();
        let mut inference = InterfaceInference::default();
        for head in [SPENSO_TAG.dot, ETS.metric] {
            let dot = FunctionBuilder::new(head).add_arg(&p).add_arg(&q).finish();
            for _ in 0..2 {
                assert!(
                    inference
                        .infer_validated(dot.as_view())
                        .unwrap()
                        .canonical()
                        .is_scalar()
                );
            }
            assert!(
                inference
                    .reusable_interfaces
                    .contains_key(dot.as_view().get_data())
            );
            let external = inference.infer_validated((&dot * &p).as_view()).unwrap();
            assert_eq!(
                external.logical_slots(),
                vec![representation.slot(PartialIndex::open(0))]
            );

            let invalid = FunctionBuilder::new(head)
                .add_arg(&p)
                .add_arg(&incompatible)
                .finish();
            assert!(inference.infer_validated(invalid.as_view()).is_err());
            assert!(
                !inference
                    .reusable_interfaces
                    .contains_key(invalid.as_view().get_data())
            );
        }
    }

    #[test]
    fn scalar_dot_inference_retains_operand_materialization_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let callback = spenso::tensor_symbol!(
            "cached_dot_materialization_callback",
            norm = move |node, _| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| {
                        Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok()
                    })
                {
                    observed.fetch_add(1, Ordering::Relaxed);
                }
            }
        );
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let p = FunctionBuilder::new(callback)
            .add_arg(representation.to_symbolic([]))
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("cached_dot_plain_operand"))
            .add_arg(representation.to_symbolic([]))
            .finish();
        let dot = FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(p)
            .add_arg(q)
            .finish();
        let mut inference = InterfaceInference::default();
        calls.store(0, Ordering::Relaxed);
        for _ in 0..2 {
            let previous = calls.load(Ordering::Relaxed);
            assert!(
                inference
                    .infer_validated(dot.as_view())
                    .unwrap()
                    .canonical()
                    .is_scalar()
            );
            assert!(calls.load(Ordering::Relaxed) > previous);
            assert!(
                !inference
                    .reusable_interfaces
                    .contains_key(dot.as_view().get_data())
            );
        }
        let previous = calls.load(Ordering::Relaxed);
        SymbolicTensor::<PartialStructure>::validate_interface(
            &dot,
            &PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        assert_eq!(calls.load(Ordering::Relaxed), previous);
    }

    #[test]
    fn result_observation_keeps_encoded_ports_without_replaying_constructor_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let callback = spenso::tensor_symbol!(
            "observed_result_materialization_changes_rank",
            norm = move |node, output| {
                observed.fetch_add(1, Ordering::Relaxed);
                if let AtomView::Fun(function) = node
                    && function.iter().last().is_some_and(|argument| {
                        Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok()
                    })
                {
                    **output = Atom::Zero;
                }
            }
        );
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let metadata = FunctionBuilder::new(symbolica::symbol!(
            "observed_result_scalar_metadata"; Scalar
        ))
        .add_arg(representation.to_symbolic([]))
        .finish();
        let structure =
            PartialStructure::from_logical_slots([representation.slot(PartialIndex::open(0))]);
        for arguments in [vec![], vec![metadata]] {
            let atom = FunctionBuilder::new(callback)
                .add_args(arguments)
                .add_arg(representation.to_symbolic([]))
                .finish();
            calls.store(0, Ordering::Relaxed);
            let inferred = InterfaceInference::default()
                .infer_validated(atom.as_view())
                .unwrap();
            assert!(inferred.canonical().is_scalar());
            assert!(calls.load(Ordering::Relaxed) > 0);

            calls.store(0, Ordering::Relaxed);
            SymbolicTensor::<PartialStructure>::validate_interface(&atom, &structure).unwrap();
            assert_eq!(calls.load(Ordering::Relaxed), 0);

            let x = Atom::var(symbolica::symbol!("observed_result_x"));
            let y = Atom::var(symbolica::symbol!("observed_result_y"));
            let value =
                SymbolicTensor::<PartialStructure>::new((&x + &y) * &atom, structure.clone())
                    .unwrap();
            let expanded = value.expression.expand();
            calls.store(0, Ordering::Relaxed);
            for finish in [
                SymbolicTensor::<PartialStructure>::with_checked_expression,
                SymbolicTensor::<PartialStructure>::with_algebra_result,
            ] {
                let result = finish(&value, expanded.clone()).unwrap();
                assert_eq!(result.expression, expanded);
                assert_eq!(result.structure, structure);
                assert_eq!(calls.load(Ordering::Relaxed), 0);
            }
        }
        // Actual rank loss is still rejected; only zero may retain a typed interface.
        assert!(
            SymbolicTensor::<PartialStructure>::validate_interface(&Atom::one(), &structure)
                .is_err()
        );
        SymbolicTensor::<PartialStructure>::validate_interface(&Atom::Zero, &structure).unwrap();
    }

    #[test]
    fn tensor_power_lowering_preserves_normalizer_calls_on_unchanged_functions() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = symbolica::symbol!(
            "power_lowering_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let atom = FunctionBuilder::new(head).add_arg(3).finish();
        calls.store(0, Ordering::Relaxed);
        assert!(
            InterfaceInference::lower_tensor_powers(atom.as_view())
                .unwrap()
                .is_none()
        );
        assert!(calls.load(Ordering::Relaxed) > 0);
    }

    #[test]
    fn recursive_merging_rejects_three_compatible_explicit_ports() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(11);
        let interfaces = vec![
            explicit_interface(representation, index),
            explicit_interface(representation, index),
            explicit_interface(representation, index),
        ];

        assert!(InterfaceInference::merge_explicit_interface_sequence(&interfaces).is_err());
    }

    #[test]
    fn recursive_merging_contracts_one_explicit_pair() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(13);
        let interfaces = vec![
            explicit_interface(representation, index),
            explicit_interface(representation, index),
        ];

        assert!(
            InterfaceInference::merge_explicit_interface_sequence(&interfaces)
                .unwrap()
                .canonical()
                .is_scalar()
        );
    }

    #[test]
    fn recursive_merging_normalizes_explicit_indices_within_one_interface() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(15);
        let pair = PartialStructure::from_logical_slots([
            representation.slot(PartialIndex::Explicit(index)),
            representation.slot(PartialIndex::Explicit(index)),
        ]);
        let triple = PartialStructure::from_logical_slots([
            representation.slot(PartialIndex::Explicit(index)),
            representation.slot(PartialIndex::Explicit(index)),
            representation.slot(PartialIndex::Explicit(index)),
        ]);

        assert!(
            InterfaceInference::merge_explicit_interface_sequence(&[pair])
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(InterfaceInference::merge_explicit_interface_sequence(&[triple]).is_err());
    }

    #[test]
    fn additive_reinference_compares_canonical_scalar_interfaces() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let diagonal = |name: &str, index| {
            let slot = representation.slot::<AbstractIndex, _>(index).to_atom();
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(&slot)
                .add_arg(&slot)
                .finish()
        };
        let left = diagonal("additive_scalar_left", AbstractIndex::Normal(21));
        let right = diagonal("additive_scalar_right", AbstractIndex::Normal(23));

        assert!(
            InterfaceInference::default()
                .infer(&(left.as_ref() + Atom::num(1)))
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(
            InterfaceInference::default()
                .infer(&(left + right))
                .unwrap()
                .canonical()
                .is_scalar()
        );
    }

    #[test]
    fn powers_reject_more_than_two_explicit_index_occurrences() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(17);
        let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("powered_explicit_index"))
            .add_arg(representation.slot::<AbstractIndex, _>(index).to_atom())
            .finish();

        assert!(
            InterfaceInference::default()
                .infer(&tensor.as_ref().pow(Atom::num(2)))
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(
            InterfaceInference::default()
                .infer(&tensor.as_ref().pow(Atom::num(3)))
                .is_err()
        );
    }

    #[test]
    fn powers_count_repeated_explicit_indices_in_the_base() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(19);
        let slot = representation.slot::<AbstractIndex, _>(index).to_atom();
        let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("powered_repeated_index"))
            .add_arg(&slot)
            .add_arg(&slot)
            .finish();

        assert!(
            InterfaceInference::default()
                .infer(&tensor.as_ref().pow(Atom::num(2)))
                .is_err()
        );
    }

    #[test]
    fn open_ports_retain_power_parity_semantics() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("powered_open_index"))
            .add_arg(representation.to_symbolic([]))
            .finish();

        assert!(
            InterfaceInference::default()
                .infer(&tensor.as_ref().pow(Atom::num(2)))
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert_eq!(
            InterfaceInference::default()
                .infer(&tensor.as_ref().pow(Atom::num(3)))
                .unwrap()
                .canonical()
                .order(),
            1
        );
    }

    #[test]
    fn leaf_slots_follow_syntactic_argument_order() {
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let atom = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(
                euclidean
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(17))
                    .to_atom(),
            )
            .add_arg(
                minkowski
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(19))
                    .to_atom(),
            )
            .finish();

        let inferred = atom
            .as_view()
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
            .unwrap();
        assert_eq!(
            inferred
                .layout()
                .canonical_to_logical(&inferred.canonical().external_structure())
                .into_iter()
                .map(|slot| slot.rep())
                .collect::<Vec<_>>(),
            vec![euclidean, minkowski]
        );
    }
    #[test]
    fn typed_power_preserves_zero_and_unresolved_interfaces() {
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let p = spenso::vector_symbol!("typed_power_boundary::p");
        let scalar = |atom| {
            SymbolicTensor::from_produced_parts(atom, PartialStructure::from_logical_slots([]))
                .unwrap()
        };
        let explicit = SymbolicTensor::<PartialStructure>::infer(
            FunctionBuilder::new(p)
                .add_arg(
                    rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(83117))
                        .to_atom(),
                )
                .finish(),
        )
        .unwrap();
        let two = scalar(Atom::num(2));
        let three = scalar(Atom::num(3));
        let squared = explicit.pow(&two).unwrap();
        assert!(squared.is_scalar());
        assert_eq!(squared.expression, explicit.expression.pow(2));
        assert!(explicit.pow(&three).is_err());
        assert!(explicit.pow(&scalar(Atom::num(-1))).is_err());
        assert!(explicit.pow(&scalar(Atom::num((1, 2)))).is_err());
        assert!(explicit.pow(&explicit).is_err());
        assert_eq!(explicit.pow(&scalar(Atom::one())).unwrap(), explicit);
        assert!(explicit.pow(&scalar(Atom::Zero)).unwrap().is_scalar());

        let zero =
            SymbolicTensor::from_produced_parts(Atom::Zero, explicit.structure.clone()).unwrap();
        let squared_zero = zero.pow(&two).unwrap();
        assert!(squared_zero.is_zero());
        assert!(squared_zero.is_scalar());
        assert!(zero.pow(&three).is_err());
        let open = SymbolicTensor::<PartialStructure>::infer(
            FunctionBuilder::new(p)
                .add_arg(rep.to_symbolic([]))
                .finish(),
        )
        .unwrap();
        for (power, rank) in [(2, 0), (3, 1)] {
            let result = open.pow(&scalar(Atom::num(power))).unwrap();
            assert_eq!(result.rank(), rank);
            SymbolicTensor::validate_atom(&result.expression).unwrap();
            SymbolicTensor::validate_interface(&result.expression, &result.structure).unwrap();
        }
        let open_zero =
            SymbolicTensor::from_produced_parts(Atom::Zero, open.structure.clone()).unwrap();
        let result = open_zero.pow(&three).unwrap();
        assert!(result.is_zero());
        assert_eq!(
            result.structure.logical_slots(),
            open.structure.logical_slots()
        );

        let x = Atom::var(symbolica::symbol!("typed_power_boundary::x"));
        let base = scalar(&x + 1);
        for exponent in [Atom::num(-1), Atom::num((1, 2)), x] {
            assert_eq!(
                base.pow(&scalar(exponent.clone())).unwrap().expression,
                base.expression.pow(&exponent)
            );
        }
    }

    #[test]
    fn typed_division_keeps_scalar_contractions_scoped() {
        crate::test_support::test_initialize();
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let slot = rep
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(83117))
            .to_atom();
        let p = FunctionBuilder::new(spenso::vector_symbol!("typed_division_p"))
            .add_arg(&slot)
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("typed_division_q"))
            .add_arg(&slot)
            .finish();
        let denominator = SymbolicTensor::<PartialStructure>::infer(&p * &q)
            .unwrap()
            .finish_parts()
            .unwrap();
        assert!(denominator.is_scalar());
        let numerator = SymbolicTensor::<PartialStructure>::infer(p).unwrap();
        let result = numerator.divide(&denominator).unwrap();
        assert_eq!(
            result.structure.logical_slots(),
            numerator.structure.logical_slots()
        );
        SymbolicTensor::validate_interface(&result.expression, &result.structure).unwrap();
        SymbolicTensor::validate_atom(&result.expression).unwrap();
        let expected =
            &numerator.expression * spenso::bracket!(denominator.expression.clone()).pow(-1);
        assert_eq!(result.expression, expected);
        assert!(denominator.divide(&numerator).is_err());
    }

    #[test]
    fn tensor_value_equality_uses_logical_ports_and_ignores_cached_facts() {
        use std::hash::{DefaultHasher, Hasher};
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let slots = [83121, 83123]
            .map(|index| rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(index))));
        let atom = FunctionBuilder::new(spenso::tensor_symbol!("tensor_value_boundary"))
            .add_args(slots.map(composition::port_atom))
            .finish();
        let value =
            SymbolicTensor::new(atom.clone(), PartialStructure::from_logical_slots(slots)).unwrap();
        let mut same = value.clone();
        same.is_metric = !value.is_metric;
        same.is_composite = !value.is_composite;
        assert!(same.rewrites_preserve_interface());
        assert!(value.same_value(&same));
        let hash = |tensor: &SymbolicTensor<PartialStructure>| {
            let mut state = DefaultHasher::new();
            tensor.hash_value(&mut state);
            state.finish()
        };
        assert_eq!(hash(&value), hash(&same));
        let reversed = SymbolicTensor::new(
            atom,
            PartialStructure::from_logical_slots(slots.into_iter().rev()),
        )
        .unwrap();
        assert!(!value.same_value(&reversed));
        let scalar_zero =
            SymbolicTensor::new(Atom::Zero, PartialStructure::from_logical_slots([])).unwrap();
        let tensor_zero = SymbolicTensor::new(Atom::Zero, value.structure.clone()).unwrap();
        assert!(scalar_zero.is_zero() && tensor_zero.is_zero());
        assert!(!scalar_zero.same_value(&tensor_zero));
    }

    #[test]
    fn finishing_inferred_carrier_retains_its_validation() {
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let source = SymbolicTensor::<PartialStructure>::infer(
            FunctionBuilder::new(spenso::tensor_symbol!("finished_proof_boundary"))
                .add_arg(
                    rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(83131))
                        .to_atom(),
                )
                .finish(),
        )
        .unwrap();
        assert!(source.proofs.is_admitted());
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        SCOPE_VALIDATIONS.with(|count| count.set(0));
        let result = source.clone().finish_parts().unwrap();
        assert!(result.same_value(&source));
        assert!(result.rewrites_preserve_interface());
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
        assert!(
            SymbolicTensor::from_declared_parts(
                result.expression,
                PartialStructure::from_logical_slots([]),
            )
            .is_err()
        );
    }
    #[test]
    fn canonical_rank_one_components_have_scalar_interfaces() {
        crate::test_support::test_initialize();
        let vector = spenso::vector_symbol!("typed_component::vector");
        let metadata = Atom::var(symbolica::symbol!("typed_component::edge"));
        for marker in [AIND_SYMBOLS.cind, AIND_SYMBOLS.find] {
            let component = FunctionBuilder::new(vector)
                .add_arg(&metadata)
                .add_arg(symbolica::function!(marker, 0))
                .finish();
            let typed = SymbolicTensor::infer(component.clone()).unwrap();
            assert!(typed.structure().canonical().is_scalar());
            assert_eq!(typed.expression(), &component);
            let factored = (Atom::one() + &metadata) * &component;
            let typed = SymbolicTensor::infer(factored.clone()).unwrap();
            assert!(typed.structure().canonical().is_scalar());
            assert_eq!(typed.expression(), &factored);
        }
    }

    #[test]
    fn delta_component_selector_is_metadata_and_contracts_with_minkowski_signs() {
        crate::test_support::test_initialize();
        let vector = spenso::vector_symbol!("typed_component::delta_momentum");
        let slot = ExtendibleReps::MINKOWSKI
            .new_rep(Dimension::Concrete(4))
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(881924))
            .to_atom();
        for component in 0..4 {
            let selector = symbolica::function!(AIND_SYMBOLS.cind, component);
            let delta = symbolica::function!(ETS.delta, &selector, &slot);
            let typed = SymbolicTensor::infer(delta.clone()).unwrap();
            assert_eq!(typed.structure().canonical().order(), 1);
            assert_eq!(typed.expression(), &delta);
            let contracted = SymbolicTensor::infer(delta * symbolica::function!(vector, &slot))
                .unwrap()
                .contract(Default::default())
                .unwrap()
                .expand_dots()
                .unwrap();
            let sign = if component == 0 { 1 } else { -1 };
            assert_eq!(
                contracted.expression(),
                &(Atom::num(sign) * symbolica::function!(vector, selector))
            );
        }
        let invalid =
            symbolica::function!(ETS.delta, symbolica::function!(AIND_SYMBOLS.cind, -1), slot);
        assert!(SymbolicTensor::infer(invalid).is_err());
    }

    #[test]
    fn nested_fixed_vectors_are_admitted_before_scalar_powers_and_scopes() {
        crate::test_support::test_initialize();
        let momentum = spenso::vector_symbol!("nested_admission::p");
        let metadata = symbolica::function!(
            symbolica::symbol!("nested_admission::provenance"; Scalar),
            6,
            -symbolica::function!(symbolica::symbol!("nested_admission::source"), 2)
        );
        let dimension = Atom::var(symbolica::symbol!("nested_admission::d"));
        let x = Atom::var(symbolica::symbol!("nested_admission::x"; Scalar));
        let y = Atom::var(symbolica::symbol!("nested_admission::y"; Scalar));
        let spectator = (&x + &y).pow(9)
            * symbolica::function!(
                symbolica::symbol!("nested_admission::opaque"; Scalar),
                (&x - &y).pow(11)
            )
            / (&x + Atom::num(3));
        for dimension in [dimension, Atom::num(4)] {
            let rep = Minkowski {}.new_rep(Dimension::try_from(dimension.as_view()).unwrap());
            let compact = rep.to_symbolic([]);
            let slot = rep
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(881925))
                .to_atom();
            let vector = symbolica::function!(momentum, &metadata, &compact);
            let generic = symbolica::function!(momentum, &metadata, &vector);
            assert_eq!(
                SymbolicTensor::infer(generic)
                    .unwrap()
                    .to_dots()
                    .unwrap()
                    .expression(),
                &spenso::dot!(&vector, &vector)
            );
            for component in 0..4 {
                let coordinate = symbolica::function!(AIND_SYMBOLS.cind, component);
                let nested = symbolica::function!(ETS.delta, &coordinate, &vector);
                let expected = spenso::dot!(
                    symbolica::function!(ETS.delta, &coordinate, &compact),
                    &vector
                );
                let admitted = SymbolicTensor::infer(nested.clone()).unwrap();
                assert!(admitted.is_scalar());
                assert_eq!(admitted.to_dots().unwrap().expression(), &expected);

                // Explicit repeated indices independently define the contraction.
                let indexed = SymbolicTensor::infer(
                    symbolica::function!(ETS.delta, &coordinate, &slot)
                        * symbolica::function!(momentum, &metadata, &slot),
                )
                .unwrap()
                .contract(Default::default())
                .unwrap()
                .to_dots()
                .unwrap();
                assert_eq!(indexed.expression(), &expected);
                for basis in 0..4 {
                    let fixed = symbolica::function!(
                        ETS.delta,
                        symbolica::function!(AIND_SYMBOLS.cind, basis),
                        &compact
                    );
                    let actual = admitted
                        .to_dots()
                        .unwrap()
                        .expression()
                        .replace(vector.clone())
                        .with(fixed.clone());
                    let nested_basis =
                        SymbolicTensor::infer(symbolica::function!(ETS.delta, &coordinate, fixed))
                            .unwrap();
                    assert_eq!(nested_basis.expression(), &actual);
                    assert_eq!(
                        actual,
                        Atom::num(if component != basis {
                            0
                        } else if component == 0 {
                            1
                        } else {
                            -1
                        })
                    );
                }
                for power in [-2, 1, 2, 3] {
                    let source = &spectator * (nested.pow(power) + &x);
                    let result = SymbolicTensor::infer(source).unwrap().to_dots().unwrap();
                    assert!(result.is_scalar());
                    assert_eq!(
                        result.expression(),
                        &(&spectator * (expected.pow(power) + &x))
                    );
                    assert_eq!(
                        SymbolicTensor::infer(result.expression().clone()).unwrap(),
                        result
                    );
                }
                let open = symbolica::function!(momentum, &metadata, &slot);
                let product = SymbolicTensor::infer(&spectator * &nested * &open)
                    .unwrap()
                    .to_dots()
                    .unwrap();
                assert_eq!(product.expression(), &(&spectator * &expected * &open));
                assert_eq!(
                    product.structure(),
                    SymbolicTensor::infer(open).unwrap().structure()
                );
            }
        }
    }

    #[test]
    fn nested_vector_admission_retains_callback_and_metadata_guards() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };

        crate::test_support::test_initialize();
        let vector = spenso::vector_symbol!("nested_admission_guard::p");
        let rep = Minkowski {}.new_rep(4);
        let compact = rep.to_symbolic([]);
        let slot = rep
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(881926))
            .to_atom();
        let coordinate = symbolica::function!(AIND_SYMBOLS.cind, 0);
        let p = symbolica::function!(vector, &compact);
        for invalid in [
            symbolica::function!(vector, &slot),
            symbolica::function!(
                spenso::tensor_symbol!("nested_admission_guard::rank_two"),
                &compact,
                &compact
            ),
            symbolica::function!(vector, &compact, &compact),
            symbolica::function!(vector, symbolica::function!(vector, &slot), &compact),
        ] {
            assert!(
                SymbolicTensor::infer(symbolica::function!(ETS.delta, &coordinate, invalid))
                    .is_err()
            );
        }

        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let callback = spenso::vector_symbol!(
            "nested_admission_guard::callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let source = symbolica::function!(callback, &p);
        calls.store(0, Ordering::Relaxed);
        assert!(SymbolicTensor::infer(source).is_err());
        assert_eq!(
            calls.load(Ordering::Relaxed),
            1,
            "only the established power-lowering constructor may replay the callback"
        );

        let observed = Arc::clone(&calls);
        let ancestor = spenso::tensor_symbol!(
            "nested_admission_guard::callback_ancestor",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let source = symbolica::function!(
            ancestor,
            symbolica::function!(ETS.delta, &coordinate, &p),
            &slot
        );
        calls.store(0, Ordering::Relaxed);
        assert!(SymbolicTensor::infer(source).is_err());
        assert_eq!(
            calls.load(Ordering::Relaxed),
            1,
            "admission must not rebuild a callback-bearing ancestor"
        );

        let observed = Arc::clone(&calls);
        let metadata = symbolica::symbol!("nested_admission_guard::metadata"; Scalar;
            norm = move |_, _| { observed.fetch_add(1, Ordering::Relaxed); }
        );
        let source = symbolica::function!(
            ETS.delta,
            &coordinate,
            symbolica::function!(vector, symbolica::function!(metadata, 7), &compact)
        );
        calls.store(0, Ordering::Relaxed);
        assert!(SymbolicTensor::infer(source).is_err());
        assert_eq!(
            calls.load(Ordering::Relaxed),
            0,
            "opaque scalar metadata must not be replayed"
        );

        let scalar = symbolica::symbol!("nested_admission_guard::scalar"; Scalar);
        for unchanged in [
            spenso::g!(&p, &p),
            spenso::dot!(&p, &p),
            symbolica::function!(vector, &slot).pow(2),
            symbolica::function!(scalar, symbolica::function!(ETS.delta, coordinate, &p)),
        ] {
            assert_eq!(
                SymbolicTensor::infer(unchanged.clone())
                    .unwrap()
                    .expression(),
                &unchanged
            );
        }
    }

    #[test]
    fn concrete_rank_one_admission_keeps_shape_metadata_and_namespace_checks() {
        crate::test_support::test_initialize();
        let vector = spenso::vector_symbol!("typed_component::shape_vector");
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let slot = rep
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(881923))
            .to_atom();
        let fixed = symbolica::function!(AIND_SYMBOLS.cind, 0);
        let make = |args: Vec<Atom>| {
            let mut leaf = FunctionBuilder::new(vector);
            for argument in args {
                leaf = leaf.add_arg(argument);
            }
            leaf.finish()
        };
        for bad in [
            make(vec![]),
            make(vec![Atom::num(0)]),
            make(vec![symbolica::function!(
                symbolica::symbol!("typed_component::cind"),
                0
            )]),
            make(vec![FunctionBuilder::new(AIND_SYMBOLS.cind).finish()]),
            make(vec![symbolica::function!(AIND_SYMBOLS.cind, -1)]),
            make(vec![symbolica::function!(
                AIND_SYMBOLS.cind,
                Atom::num(1) / Atom::num(2)
            )]),
            make(vec![symbolica::function!(AIND_SYMBOLS.cind, 0, 1)]),
            make(vec![fixed.clone(), Atom::num(1)]),
            make(vec![fixed.clone(), fixed.clone()]),
            make(vec![slot.clone(), fixed.clone()]),
            make(vec![fixed.clone(), slot.clone()]),
            make(vec![
                symbolica::function!(
                    spenso::tensor_symbol!("typed_component::open_metadata"),
                    &slot
                ),
                fixed.clone(),
            ]),
        ] {
            assert!(SymbolicTensor::infer(bad.clone()).is_err(), "{bad}");
        }
        let explicit = SymbolicTensor::infer(make(vec![slot])).unwrap();
        assert_eq!(explicit.structure().canonical().order(), 1);
        let compact = SymbolicTensor::infer(make(vec![rep.to_symbolic([])])).unwrap();
        assert_eq!(compact.structure().canonical().order(), 1);
    }

    #[test]
    fn concrete_component_observation_does_not_replay_or_hide_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let vector = spenso::vector_symbol!(
            "typed_component::callback",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let component = symbolica::function!(vector, 7, symbolica::function!(AIND_SYMBOLS.cind, 2));
        calls.store(0, Ordering::Relaxed);
        let mut inference = InterfaceInference {
            leaf_inference: LeafInference::Observe,
            ..Default::default()
        };
        assert!(
            inference
                .infer_validated(component.as_view())
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(!inference.normalization_is_intrinsic_cached(component.as_view()));
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "normalization callback contract violated in with_renamed_symbols")]
    fn symbol_renaming_checks_callback_contract() {
        test_support::test_initialize();
        let source_index = symbol!("rename_callback_source_violation");
        let target_index = symbol!("rename_callback_target_violation");
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = tensor_symbol!(
            "rename_callback_tensor_violation",
            norm = move |node, out| {
                observed.fetch_add(1, Ordering::Relaxed);
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(slot)) = function.iter().next()
                    && slot.iter().nth(1) == Some(Atom::var(target_index).as_view())
                {
                    **out = Atom::one();
                }
            }
        );
        let source = SymbolicTensor::infer(
            FunctionBuilder::new(head)
                .add_arg(
                    ExtendibleReps::MINKOWSKI
                        .new_rep(4)
                        .to_symbolic([Atom::var(source_index)]),
                )
                .finish(),
        )
        .unwrap();
        calls.store(0, Ordering::Relaxed);
        source
            .rename_symbols(&HashMap::from([(source_index, target_index)]))
            .unwrap();
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(
        expected = "normalization callback contract violated in map_representation_dimensions"
    )]
    fn dimension_mapping_checks_callback_contract() {
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let head = tensor_symbol!(
            "dimension_result_normalizer_violation",
            norm = move |node, output| {
                observed.lock().unwrap().push(node.to_owned());
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(slot)) = function.iter().last()
                    && slot.iter().next() == Some(Atom::num(6).as_view())
                {
                    let mode = function.iter().next().unwrap();
                    if mode == Atom::Zero.as_view() || mode == Atom::one().as_view() {
                        **output = mode.to_owned();
                    }
                }
            }
        );
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        for index in [
            PartialIndex::open(0),
            PartialIndex::Explicit(AbstractIndex::Normal(73321)),
        ] {
            for mode in [0, 1, 2] {
                let port = representation.slot(index);
                let expression = FunctionBuilder::new(head)
                    .add_arg(mode)
                    .add_arg(composition::port_atom(port))
                    .finish();
                let value = SymbolicTensor::<PartialStructure>::new(
                    expression,
                    PartialStructure::from_logical_slots([port]),
                )
                .unwrap();
                calls.lock().unwrap().clear();
                let unchanged = value
                    .expression
                    .with_lorentz_dimension(Atom::num(4).as_view());
                let unchanged_calls = mem::take(&mut *calls.lock().unwrap());
                let mut mapped_pairs = Vec::new();
                let actual = value
                    .map_representation_dimensions(|representation, dimension| {
                        mapped_pairs.push((representation, dimension.to_owned()));
                        dimension.to_owned()
                    })
                    .unwrap();
                assert_eq!(actual.expression, unchanged);
                assert_eq!(actual.structure, value.structure);
                assert_eq!(*calls.lock().unwrap(), unchanged_calls);
                assert_eq!(mapped_pairs, vec![(representation.rep, Atom::num(4))]);

                calls.lock().unwrap().clear();
                let expected = value
                    .expression
                    .with_lorentz_dimension(Atom::num(6).as_view());
                let expected_calls = mem::take(&mut *calls.lock().unwrap());
                assert!(!expected_calls.is_empty());
                let result = value.with_lorentz_dimension(Dimension::Concrete(6));
                assert_eq!(*calls.lock().unwrap(), expected_calls);
                if mode == 1 {
                    assert_eq!(expected, Atom::one());
                }
                let result = result.unwrap();
                assert_eq!(result.expression, expected);
                assert_eq!(
                    result.structure.logical_slots(),
                    vec![
                        ExtendibleReps::MINKOWSKI
                            .new_rep(Dimension::Concrete(6))
                            .slot(index)
                    ]
                );
            }
        }

        let tensor = FunctionBuilder::new(tensor_symbol!("dimension_promoted_tensor_violation"))
            .add_arg(ExtendibleReps::MINKOWSKI.new_rep(6).to_symbolic([]))
            .finish();
        let promoted_tensor = tensor.clone();
        let observed = Arc::clone(&calls);
        let scalar = symbol!(
            "dimension_scalar_callback_violation"; Scalar;
            norm = move |node, output| {
                observed.lock().unwrap().push(node.to_owned());
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(slot)) = function.iter().next()
                    && slot.iter().next() == Some(Atom::num(6).as_view())
                {
                    **output = promoted_tensor.clone();
                }
            }
        );
        let expression = FunctionBuilder::new(scalar)
            .add_arg(representation.to_symbolic([]))
            .finish();
        assert!(!InterfaceInference::has_structured_syntax(
            expression.as_view()
        ));
        let value = SymbolicTensor::<PartialStructure>::new(
            expression,
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        calls.lock().unwrap().clear();
        let expected = value
            .expression
            .with_lorentz_dimension(Atom::num(6).as_view());
        assert_eq!(expected, tensor);
        let expected_calls = mem::take(&mut *calls.lock().unwrap());
        let mapped = value
            .with_lorentz_dimension(Dimension::Concrete(6))
            .unwrap();
        assert_eq!(mapped.expression, tensor);
        assert!(mapped.is_scalar());
        assert_eq!(*calls.lock().unwrap(), expected_calls);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "normalization callback contract violated in wrap_indices")]
    fn scoping_checks_callback_contract() {
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let scope = symbol!("shared_scope_test_violation");
        let index = AbstractIndex::Normal(73201);
        let slot = rep.slot(PartialIndex::Explicit(index));
        let target = rep.slot::<AbstractIndex, _>(index.scoped(scope)).to_atom();
        let callback = tensor_symbol!(
            "scope_callback_strips_ports_violation",
            norm = move |node, output| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|arg| arg == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let value = SymbolicTensor::<PartialStructure>::new(
            FunctionBuilder::new(callback)
                .add_arg(composition::port_atom(slot))
                .finish(),
            PartialStructure::from_logical_slots([slot]),
        )
        .unwrap();
        value.wrap_indices(scope, false).unwrap();
    }

    fn cooking_callback_case(violate: bool) -> InferenceResult<SymbolicTensor<PartialStructure>> {
        let name = if violate {
            "cooking_contract_bad"
        } else {
            "cooking_contract_good"
        };
        let index_head = symbol!("cooking_contract_index", tags = [&SPENSO_TAG.index]);
        let head = tensor_symbol!(
            &name,
            norm = move |node, out| {
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(slot)) = function.iter().next()
                    && let Some(AtomView::Fun(index)) = slot.iter().nth(1)
                    && index.get_symbol() == index_head
                    && index.get_nargs() == 1
                {
                    **out = if violate {
                        Atom::one()
                    } else {
                        Atom::num(3)
                            * FunctionBuilder::new(tensor_symbol!("cooking_contract_result"))
                                .add_args(function.iter())
                                .finish()
                    };
                }
            }
        );
        // Both spellings denote the same supported native index. Cooking
        // canonicalizes local zero from head(owner, 0) to head(owner), rebuilding
        // the callback-bearing parent without changing any logical index.
        let alias = FunctionBuilder::new(index_head)
            .add_arg(73901)
            .add_arg(0)
            .finish();
        let slot = ExtendibleReps::EUCLIDEAN.new_rep(3).to_symbolic([alias]);
        let source = SymbolicTensor::infer(FunctionBuilder::new(head).add_arg(slot).finish())?;
        source.with_cooked_indices(&CookSettings::indices())
    }

    #[test]
    fn cooking_preserves_compliant_callback_interface() {
        let result = cooking_callback_case(false).unwrap();
        assert_eq!(
            SymbolicTensor::infer(result.expression.clone())
                .unwrap()
                .structure,
            result.structure
        );
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "normalization callback contract violated in with_cooked_indices")]
    fn cooking_checks_callback_contract() {
        cooking_callback_case(true).unwrap();
    }
    #[test]
    fn wrapping_preserves_compliant_callback_interface() {
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let label = symbol!(
            "wrapping_contract_index",
            tags = [&SPENSO_TAG.index],
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let source_head = tensor_symbol!("wrapping_contract_source");
        let head = tensor_symbol!(
            "wrapping_contract_good",
            norm = |node, out| {
                if let AtomView::Fun(function) = node {
                    **out = Atom::num(3) * function.iter().next().unwrap();
                }
            }
        );
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(3);
        for index in [
            AbstractIndex::Normal(73931),
            AbstractIndex::Open {
                owner: 73931,
                axis: 0,
            },
            AbstractIndex::Named(label.into(), 73931, 0),
        ] {
            let source = SymbolicTensor::infer(
                FunctionBuilder::new(source_head)
                    .add_arg(rep.to_symbolic([Atom::from(index)]))
                    .finish(),
            )
            .unwrap();
            let zero = source.with_checked_expression(Atom::Zero).unwrap();
            let raw = SymbolicTensor {
                structure: OrderedStructure::new(vec![rep.slot::<AbstractIndex, _>(index)])
                    .into_canonical(),
                expression: source.expression.clone(),
                is_metric: false,
                is_composite: false,
                work: Default::default(),
                proofs: TensorProofs::default(),
            };
            let raw_zero = SymbolicTensor {
                structure: raw.structure.clone(),
                expression: Atom::Zero,
                is_metric: false,
                is_composite: true,
                work: Default::default(),
                proofs: TensorProofs::default(),
            };
            let indexed = (!matches!(index, AbstractIndex::Open { .. }))
                .then(|| source.clone().try_into_indexed().unwrap());
            calls.store(0, Ordering::Relaxed);
            let wrapped_zero = Wrap {}.apply(&head, zero).unwrap();
            assert!(wrapped_zero.is_zero());
            assert_eq!(wrapped_zero.structure, source.structure);
            let raw_zero = Wrap {}.apply(&head, raw_zero).unwrap();
            assert!(raw_zero.expression.is_zero());
            assert_eq!(raw_zero.structure.order(), 1);
            let raw = Wrap {}.apply(&head, raw).unwrap();
            let indexed = indexed.map(|source| Wrap {}.apply(&head, source).unwrap());
            let wrapped = Wrap {}.apply(&head, source).unwrap();
            assert_eq!(raw.expression, wrapped.expression);
            assert_eq!(
                calls.load(Ordering::Relaxed),
                0,
                "wrapping checks must not replay index normalizers"
            );
            let inferred = SymbolicTensor::infer(wrapped.expression.clone()).unwrap();
            assert_eq!(inferred.structure, wrapped.structure);
            if let Some(indexed) = indexed {
                assert_eq!(
                    inferred.try_into_indexed().unwrap().structure,
                    indexed.structure
                );
            }
        }
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "normalization callback contract violated in Wrap::apply")]
    fn wrapping_checks_callback_contract() {
        let source = SymbolicTensor::infer(
            FunctionBuilder::new(tensor_symbol!("wrapping_contract_bad_source"))
                .add_arg(
                    ExtendibleReps::EUCLIDEAN
                        .new_rep(3)
                        .to_symbolic([Atom::num(73933)]),
                )
                .finish(),
        )
        .unwrap();
        let head = tensor_symbol!(
            "wrapping_contract_bad",
            norm = |_, out| {
                **out = Atom::one();
            }
        );
        let zero = source.with_checked_expression(Atom::Zero).unwrap();
        Wrap {}.apply(&head, zero).unwrap();
    }
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "tensor producer contract violated in from_produced_parts")]
    fn produced_parts_checks_callback_interface() {
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(3);
        let head = tensor_symbol!("produced_parts_contract", norm = |_, _| {});
        let expression = FunctionBuilder::new(head)
            .add_arg(rep.to_symbolic([Atom::num(73921)]))
            .finish();
        SymbolicTensor::from_produced_parts(
            expression,
            PartialStructure::from_logical_slots([
                rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(73923)))
            ]),
        )
        .unwrap();
    }
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "normalization callback contract violated in Wrap::apply")]
    fn wrapping_checks_encoded_free_index_identity() {
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(3);
        let index = AbstractIndex::Open {
            owner: 73941,
            axis: 0,
        };
        let source = SymbolicTensor::from_produced_parts(
            FunctionBuilder::new(tensor_symbol!("wrapping_encoded_source"))
                .add_arg(rep.to_symbolic([Atom::from(index)]))
                .finish(),
            PartialStructure::from_logical_slots([rep.slot(PartialIndex::Explicit(index))]),
        )
        .unwrap();
        let target = FunctionBuilder::new(tensor_symbol!("wrapping_encoded_target"))
            .add_arg(rep.to_symbolic([Atom::from(AbstractIndex::Open {
                owner: 73943,
                axis: 0,
            })]))
            .finish();
        let head = tensor_symbol!(
            "wrapping_encoded_bad",
            norm = move |_, out| {
                **out = target.clone();
            }
        );
        Wrap {}.apply(&head, source).unwrap();
    }
}
