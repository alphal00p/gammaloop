//! Shared scope and explicit-index observation. This is a syntax-only pass:
//! no normalization callback, materialization, or dummy allocation is permitted.
#[cfg(all(test, debug_assertions))]
use super::DebugInference;
#[cfg(test)]
use super::tests::{
    DEBUG_OCCURRENCE_FUNCTION_VISITS, DEBUG_SCOPE_VALIDATIONS, OCCURRENCE_FUNCTION_VISITS,
    SCOPE_VALIDATIONS,
};
use super::{InferenceResult, InterfaceInference, LeafInference, TensorInferenceError};
use crate::color::CS;
use crate::dirac::AGS;
use crate::tensor::TensorProofs;
use crate::tensor::composition;
use crate::tensor::composition::ExplicitIndexOccurrences;
use ahash::AHashMap;
use spenso::network::library::symbolic::ETS;
use spenso::network::parsing::{AtomStructureExt, ChainNestingError, StrictTensorFilter};
use spenso::network::tags::SPENSO_TAG;
use spenso::shadowing::{ANTISYM, CYCLIC, SYM};
use spenso::structure::TensorStructure;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::representation::{LibraryRep, RepName, Representation};
use spenso::structure::slot::{IsAbstractSlot, SlotMatch, SlotMatcher};
#[cfg(test)]
use std::cell::Cell;
use std::cell::OnceCell;
use std::mem;
#[cfg(test)]
use symbolica::atom::AtomCore;
use symbolica::atom::representation::FunView;
use symbolica::atom::{Atom, AtomView, Symbol};
use symbolica::domains::rational::Rational;

#[cfg(test)]
thread_local! {
    static SCOPED_INTERFACE_INFERENCES: Cell<usize> = const { Cell::new(0) };
}

#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub(in crate::tensor) struct ObservationScope {
    pub(super) owner: Option<Symbol>,
    pub(super) inside_factor: bool,
    pub(super) validate: bool,
    pub(super) count_indices: bool,
    pub(super) admit_direct_ports: bool,
    pub(super) pair_representation: Option<LibraryRep>,
}

type Cache = AHashMap<ObservationScope, AHashMap<Vec<u8>, ExplicitIndexOccurrences>>;

/// A producer owns its draft. Parser readers share the admitted source and
/// isolate any missing declarations before writing them, including sum terms
/// whose declared order can differ from their literal storage order.
#[derive(Clone, Copy)]
struct OccurrenceFacts<'a> {
    source: &'a TensorProofs,
    parser: Option<(AtomView<'a>, &'a OnceCell<TensorProofs>)>,
}

impl<'a> OccurrenceFacts<'a> {
    fn current(self) -> &'a TensorProofs {
        self.parser
            .and_then(|(_, local)| local.get())
            .unwrap_or(self.source)
    }

    fn writable(self) -> &'a TensorProofs {
        let Some((value, local)) = self.parser else {
            return self.source;
        };
        local.get_or_init(|| {
            let mut facts = self.source.clone();
            facts.publish_interfaces(value, self.source.retained_interface(value), []);
            facts
        })
    }

    fn interface(self, value: AtomView<'_>) -> Option<PartialStructure> {
        if self.parser.is_some() {
            self.current().inspect_interface(value)
        } else {
            self.source.admitted_interface(value)
        }
    }
}

impl InterfaceInference {
    /// Whether a function argument exports explicit tensor-index occurrences.
    ///
    /// This is the syntax-only boundary shared by construction observation and
    /// targeted index relabelling. Metadata stays opaque as a whole,
    /// including arithmetic whose tensor syntax is confined to a scalar exponent.
    pub(in crate::tensor) fn explicit_index_argument(
        parent: Symbol,
        argument: AtomView<'_>,
        slots: &mut SlotMatcher,
    ) -> bool {
        (parent.has_tag(&SPENSO_TAG.tensor) || composition::is_composite_head(parent, &SPENSO_TAG))
            && (slots.parse::<LibraryRep, AbstractIndex>(argument).is_ok()
                || Self::nested_index_argument(argument))
    }

    fn nested_index_argument(argument: AtomView<'_>) -> bool {
        Representation::<LibraryRep>::try_from(argument).is_ok()
            || argument.is_tensorial(StrictTensorFilter::Tagged)
            || matches!(argument, AtomView::Fun(nested)
                if nested.get_symbol() == *SYM
                    || nested.get_symbol() == *ANTISYM
                    || nested.get_symbol() == *CYCLIC)
    }

    /// Detect explicit index pairs in any additive branch without distributing it.
    ///
    /// The predicate suppresses ports of matching tensor leaves, without changing
    /// their payloads or recursing into scalar metadata. It must inspect borrowed
    /// syntax only and return the same decision for an unchanged leaf. Composite
    /// tensor scopes cannot be excluded. Original multiplicities and placeholder
    /// scopes are validated before exclusions; excluded leaves cannot make an
    /// invalid input valid. Closed powers and inverses keep their internal scope.
    pub fn has_explicit_index_pairs(
        value: AtomView<'_>,
        representation: LibraryRep,
        mut exclude_tensor_leaf: impl FnMut(AtomView<'_>) -> bool,
    ) -> InferenceResult<bool> {
        // Exclusions are a pair-discovery policy, not permission to hide an
        // invalid original scope. Reuse the same syntax-only observer first
        // without exclusions. Pair-query mode deliberately permits unlike
        // additive interfaces while retaining multiplicity/scope validation.
        Self::observe_indices(
            value,
            &mut SlotMatcher::default(),
            &mut AHashMap::new(),
            ObservationScope {
                validate: true,
                count_indices: true,
                pair_representation: Some(representation),
                ..ObservationScope::default()
            },
            &mut |_| false,
        )?
        .validate()?;
        // A fresh cache is required: exclusions change the counted ports.
        let occurrences = Self::observe_indices(
            value,
            &mut SlotMatcher::default(),
            &mut AHashMap::new(),
            ObservationScope {
                validate: true,
                count_indices: true,
                pair_representation: Some(representation),
                ..ObservationScope::default()
            },
            &mut exclude_tensor_leaf,
        )?;
        occurrences.validate()?;
        Ok(occurrences.selected_pair)
    }

    pub(in crate::tensor) fn index_occurrences(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        cache: &mut Cache,
    ) -> ExplicitIndexOccurrences {
        Self::observe_indices(
            value,
            slots,
            cache,
            ObservationScope {
                count_indices: true,
                ..ObservationScope::default()
            },
            &mut |_| false,
        )
        .expect("index-only observation cannot reject a placeholder scope")
    }

    pub(in crate::tensor) fn scoped_occurrences(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        #[cfg(test)]
        {
            #[cfg(debug_assertions)]
            let debug = DebugInference::is_active();
            #[cfg(not(debug_assertions))]
            let debug = false;
            if debug {
                DEBUG_SCOPE_VALIDATIONS.with(|count| count.set(count.get() + 1));
            } else {
                SCOPE_VALIDATIONS.with(|count| count.set(count.get() + 1));
            }
        }
        Self::observe_indices(
            value,
            slots,
            &mut AHashMap::new(),
            ObservationScope {
                validate: true,
                count_indices: true,
                ..ObservationScope::default()
            },
            &mut |_| false,
        )
    }

    /// An admitted metric whose arguments are tensors consumes one axis from
    /// each. A direct slot, representation or fixed component belongs to the
    /// ordinary metric constructor. Only the two argument heads are inspected.
    pub(in crate::tensor) fn is_compound_metric(
        function: FunView<'_>,
        slots: &mut SlotMatcher,
    ) -> bool {
        function.get_symbol() == ETS.metric
            && function.get_nargs() == 2
            && function.iter().all(|argument| {
                !matches!(argument, AtomView::Num(_) | AtomView::Var(_))
                    && matches!(slots.classify(argument), SlotMatch::Other)
                    && slots.concrete_component(argument).is_none()
            })
    }

    /// Constructor signatures establish which arguments are raw graph ports.
    /// This also supplies the ports for a newly emitted factor in a Dirac word.
    pub(in crate::tensor) fn constructor_ports(
        function: FunView<'_>,
        slots: &mut SlotMatcher,
        admit_direct_ports: bool,
    ) -> Option<Vec<Atom>> {
        let head = function.get_symbol();
        if head.has_tag(&SPENSO_TAG.rank1) {
            let port = slots.vector_argument(function)?;
            // A selected concrete component has no external port. Momentum
            // labels and other scalar metadata are never parser incidences.
            if slots
                .concrete_component(port)
                .is_some_and(|component| component.is_ok())
            {
                return Some(Vec::new());
            }
            return (Self::is_chain_placeholder(port)
                || slots.parse::<LibraryRep, AbstractIndex>(port).is_ok()
                || slots.parse_representation::<LibraryRep>(port).is_ok())
            .then(|| vec![port.to_owned()]);
        }
        if Self::is_compound_metric(function, slots) {
            // A metric of compound operands consumes just one axis from each.
            // Their spectators belong to the composed boundary, not an empty leaf.
            return None;
        }
        if head == CS.d {
            // The symmetric invariant d_R(R, a, b, ...) declares its group
            // representation first. The three-index d(a, b, c) has no prefix.
            return Some(
                function
                    .iter()
                    .skip(usize::from(function.get_nargs() >= 4))
                    .map(|port| port.to_owned())
                    .collect(),
            );
        }
        (head == AGS.gamma0
            || Self::is_builtin_tensor(head)
            || (admit_direct_ports && composition::is_tensor_leaf_head(head)))
        .then(|| {
            function
                .iter()
                .filter(|&port| {
                    Self::is_chain_placeholder(port)
                        || slots.parse::<LibraryRep, AbstractIndex>(port).is_ok()
                        || slots.parse_representation::<LibraryRep>(port).is_ok()
                })
                .map(|port| port.to_owned())
                .collect()
        })
    }

    /// Account for emitted indices using the producer's retained interfaces.
    /// Missing scoped facts are an internal producer error; this path never
    /// falls back to observing the free interface of an emitted expression.
    pub(in crate::tensor) fn produced_occurrences(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        proofs: &TensorProofs,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        Self::constructor_occurrences(value, slots, proofs, false)
    }

    /// At the public parts boundary the raw expression is supplied for the
    /// first time, so record the written ports of tagged generic constructors.
    pub(in crate::tensor) fn constructor_occurrences(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        proofs: &TensorProofs,
        admit_direct_ports: bool,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        Self::occurrences_with_facts(
            value,
            slots,
            OccurrenceFacts {
                source: proofs,
                parser: None,
            },
            admit_direct_ports,
        )
    }

    /// An admitted parser may need a provisional boundary for a newly scoped
    /// leaf. Return any private declarations to the caller publishing that leaf.
    pub(in crate::tensor) fn parser_occurrences(
        value: AtomView<'_>,
        source: &TensorProofs,
        boundary: Option<&PartialStructure>,
    ) -> InferenceResult<(ExplicitIndexOccurrences, Option<TensorProofs>)> {
        let local = OnceCell::new();
        let facts = OccurrenceFacts {
            source,
            parser: Some((value, &local)),
        };
        if let Some(boundary) = boundary
            && source.retained_interface(value).is_none()
        {
            facts.writable().supply_interface(value, boundary);
        }
        let history =
            Self::occurrences_with_facts(value, &mut SlotMatcher::default(), facts, false)?;
        Ok((history, local.into_inner()))
    }

    fn occurrences_with_facts(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        facts: OccurrenceFacts<'_>,
        admit_direct_ports: bool,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        // The per-call cache can reuse children. Its root is returned directly.
        Self::observe_indices_uncached(
            value,
            slots,
            &mut AHashMap::new(),
            ObservationScope {
                validate: true,
                admit_direct_ports,
                count_indices: true,
                ..ObservationScope::default()
            },
            &mut |_| false,
            Some(facts),
        )
    }

    /// Keep only the boundary ports of one independently scoped expression.
    /// Observing its interface cannot normalize callbacks or allocate indices.
    fn scope_explicit_occurrences(
        value: AtomView<'_>,
        occurrences: &mut ExplicitIndexOccurrences,
        slots: &mut SlotMatcher,
        proofs: Option<OccurrenceFacts<'_>>,
    ) -> InferenceResult<PartialStructure> {
        let interface = if let Some(proofs) = proofs {
            proofs.interface(value).ok_or_else(|| {
                TensorInferenceError::invalid(format!(
                    "internal producer omitted the scoped tensor interface for {value}"
                ))
            })?
        } else {
            #[cfg(test)]
            {
                #[cfg(debug_assertions)]
                let debug = DebugInference::is_active();
                #[cfg(not(debug_assertions))]
                let debug = false;
                if !debug {
                    SCOPED_INTERFACE_INFERENCES.with(|count| count.set(count.get() + 1));
                }
            }
            let mut inference = InterfaceInference {
                slots: mem::take(slots),
                leaf_inference: LeafInference::Observe,
                ..InterfaceInference::default()
            };
            let interface = inference
                .infer_view(value)
                .and_then(|interface| Self::merge_explicit_interface_sequence(&[interface]));
            *slots = inference.slots;
            interface?
        };
        let boundary = interface.logical_slots();
        occurrences.retain(|slot| {
            boundary.iter().any(|port| {
                port.rep().base() == slot.rep().base()
                    && port.aind == PartialIndex::Explicit(slot.aind())
            })
        });
        Ok(interface)
    }

    pub(super) fn observe_indices(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        cache: &mut Cache,
        scope: ObservationScope,
        exclude_tensor_leaf: &mut impl FnMut(AtomView<'_>) -> bool,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        Self::observe_indices_with_interfaces(value, slots, cache, scope, exclude_tensor_leaf, None)
    }

    fn observe_indices_with_interfaces(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        cache: &mut Cache,
        scope: ObservationScope,
        exclude_tensor_leaf: &mut impl FnMut(AtomView<'_>) -> bool,
        proofs: Option<OccurrenceFacts<'_>>,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        let reusable = matches!(value, AtomView::Fun(_))
            && !value.needs_normalization()
            && value.get_byte_size() <= Self::CACHE_KEY_BYTES;
        if reusable
            && let Some(summary) = cache
                .get(&scope)
                .and_then(|values| values.get(value.get_data()))
        {
            return Ok(summary.clone());
        }
        let result = Self::observe_indices_uncached(
            value,
            slots,
            cache,
            scope,
            exclude_tensor_leaf,
            proofs,
        )?;
        if reusable
            && cache.values().map(|values| values.len()).sum::<usize>() < Self::CACHE_ENTRIES
        {
            // Scope is part of the key: a placeholder proven inside a chain is
            // not legal when the same literal subtree occurs outside it.
            cache
                .entry(scope)
                .or_default()
                .insert(value.get_data().to_vec(), result.clone());
        }
        Ok(result)
    }

    fn observe_indices_uncached(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        cache: &mut Cache,
        mut scope: ObservationScope,
        exclude_tensor_leaf: &mut impl FnMut(AtomView<'_>) -> bool,
        proofs: Option<OccurrenceFacts<'_>>,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        #[cfg(test)]
        if scope.count_indices && matches!(value, AtomView::Fun(_)) {
            #[cfg(debug_assertions)]
            let debug = DebugInference::is_active();
            #[cfg(not(debug_assertions))]
            let debug = false;
            if debug {
                DEBUG_OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(count.get() + 1));
            } else {
                OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(count.get() + 1));
            }
        }
        let mut result = ExplicitIndexOccurrences::default();
        if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(value) {
            if scope.count_indices {
                result.add(slot, 1);
            }
            // Admission owns the complete index identity. Names such as `in`
            // or `out` inside it are data, not chain wiring placeholders.
            return Ok(result);
        }
        if slots.parse_representation::<LibraryRep>(value).is_ok() {
            // A compact representation similarly owns its dimension metadata.
            return Ok(result);
        }
        if !scope.count_indices && !scope.validate {
            return Ok(result);
        }
        match value {
            AtomView::Var(_) if scope.validate && Self::is_chain_placeholder(value) => {
                return Err(TensorInferenceError::invalid(
                    "Spenso chain placeholders are only valid inside chain or trace tensor factors",
                ));
            }
            AtomView::Add(sum) => {
                // A producer declaring a sum's boundary also declares that
                // boundary for each alternative. Transport that contract to
                // fresh terms; powers cannot make the same deduction for bases.
                if let Some(proofs) = proofs
                    && let Some(interface) = proofs.current().retained_interface(value)
                {
                    for term in sum.iter() {
                        if proofs.current().retained_interface(term).is_none() {
                            proofs.writable().supply_interface(term, &interface);
                        }
                    }
                }
                for term in sum.iter() {
                    let mut branch = Self::observe_indices_with_interfaces(
                        term,
                        slots,
                        cache,
                        scope,
                        exclude_tensor_leaf,
                        proofs,
                    )?;
                    // Each alternative owns its internal contractions. Only
                    // surviving ports can meet factors outside this sum. Keep
                    // malformed local counts visible to the existing validator.
                    // Pair queries allow unlike branch interfaces: scope only
                    // proven interfaces there, without adding typed admission.
                    if !branch.is_empty()
                        && branch.validate().is_ok()
                        && let Err(error) =
                            Self::scope_explicit_occurrences(term, &mut branch, slots, proofs)
                        && scope.validate
                        && scope.pair_representation.is_none()
                    {
                        return Err(error);
                    }
                    result.branch_max(&branch);
                }
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    result.append(&Self::observe_indices_with_interfaces(
                        factor,
                        slots,
                        cache,
                        scope,
                        exclude_tensor_leaf,
                        proofs,
                    )?);
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                // A closed base owns its dummy indices for every power. An
                // open base repeats its explicit ports, so those counts remain
                // visible to the surrounding product. Observe never normalizes
                // leaves or allocates dummy indices.
                for (position, term) in [base, exponent].into_iter().enumerate() {
                    let child = if position == 1 {
                        ObservationScope {
                            owner: None,
                            inside_factor: false,
                            ..scope
                        }
                    } else {
                        scope
                    };
                    let mut branch = Self::observe_indices_with_interfaces(
                        term,
                        slots,
                        cache,
                        child,
                        exclude_tensor_leaf,
                        proofs,
                    )?;
                    if position == 0 && !branch.is_empty() && branch.validate().is_ok() {
                        result.selected_pair |= branch.selected_pair;
                        let interface = match Self::scope_explicit_occurrences(
                            base,
                            &mut branch,
                            slots,
                            proofs,
                        ) {
                            Ok(interface) => Some(interface),
                            Err(error) if scope.validate => return Err(error),
                            Err(_) => None,
                        };
                        if interface
                            .as_ref()
                            .is_some_and(|interface| interface.canonical().is_scalar())
                        {
                            continue;
                        }
                        if scope.validate
                            && let Some(interface) = &interface
                        {
                            // Symbolica can spell two typed, same-oriented factors
                            // as a power. Only that exact admitted product can
                            // certify their two uncontracted copies.
                            if !proofs.is_some_and(|proofs| {
                                Self::admitted_oriented_square(
                                    proofs.current(),
                                    value,
                                    interface,
                                    exponent,
                                )
                                .is_some()
                            }) {
                                Self::power_interface(interface, exponent, || value.to_string())?;
                            }
                        }
                        if let Ok(exponent) = Rational::try_from(exponent)
                            && exponent.denominator() == 1
                        {
                            let repetitions = exponent.numerator().abs();
                            let repetitions = usize::try_from(repetitions).unwrap_or(3).min(3);
                            branch.repeat(repetitions);
                        }
                    }
                    if position == 1 && branch.validate().is_ok() {
                        // The exponent is an independent scalar scope. Keep
                        // invalid local multiplicities visible to validation,
                        // but never expose its closed contractions outside it.
                        result.selected_pair |= branch.selected_pair;
                        if !Self::has_structured_syntax(exponent) {
                            continue;
                        }
                        let interface = match Self::scope_explicit_occurrences(
                            exponent,
                            &mut branch,
                            slots,
                            proofs,
                        ) {
                            Ok(interface) => Some(interface),
                            Err(error) if scope.validate => return Err(error),
                            Err(_) => None,
                        };
                        if interface
                            .as_ref()
                            .is_some_and(|interface| interface.canonical().is_scalar())
                        {
                            continue;
                        }
                        if scope.validate && interface.is_some() {
                            return Err(TensorInferenceError::invalid(
                                "a tensor exponent must be scalar",
                            ));
                        }
                    }
                    result.branch_max(&branch);
                }
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                let tags = &*SPENSO_TAG;
                if let Some(proofs) = proofs {
                    // New known constructors publish their raw declaration in
                    // this existing producer walk. A generic tensor's root
                    // boundary cannot establish hidden or ordered raw ports.
                    let recorded = proofs.current().leaf_interfaces.contains(value);
                    if !recorded
                        && let Some(ports) =
                            Self::constructor_ports(function, slots, scope.admit_direct_ports)
                    {
                        // Only publish the constructor boundary when it is
                        // not already the producer's authoritative root.
                        let proofs = proofs.writable();
                        let retained = proofs.retained_interface(value);
                        proofs.retain_leaf_ports(value, ports)?;
                        if let Some(retained) = retained {
                            proofs.retain_boundary(value, &retained);
                        }
                    }
                }
                scope.validate &= !tags.is_scalar_metadata(symbol);
                if scope.validate {
                    scope.owner = ChainNestingError::enter(scope.owner, symbol, tags)
                        .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
                }
                let tensor_leaf = symbol.has_tag(&tags.tensor);
                let count_arguments = scope.count_indices
                    && (tensor_leaf || composition::is_composite_head(symbol, tags))
                    && !(composition::is_tensor_leaf_head(symbol) && exclude_tensor_leaf(value));
                let mut inputs = 0;
                let mut outputs = 0;
                for (position, argument) in function.iter().enumerate() {
                    let mut child = scope;
                    if symbol == tags.chain {
                        child.inside_factor = position >= 2;
                    } else if symbol == tags.trace {
                        child.inside_factor = position >= 1;
                    } else if scope.validate && Self::is_chain_placeholder(argument) {
                        if !(scope.inside_factor && tensor_leaf) {
                            return Err(TensorInferenceError::invalid(
                                "Spenso chain placeholders are only valid as tensor ports inside chain or trace factors",
                            ));
                        }
                        let AtomView::Var(variable) = argument else {
                            unreachable!()
                        };
                        if variable.get_symbol() == tags.chain_in {
                            inputs += 1;
                        } else {
                            outputs += 1;
                        }
                        continue;
                    }
                    match argument {
                        AtomView::Num(_) => continue,
                        AtomView::Var(variable)
                            if !child.validate
                                || ![tags.chain_in, tags.chain_out]
                                    .contains(&variable.get_symbol()) =>
                        {
                            continue;
                        }
                        _ => {}
                    }
                    if count_arguments
                        && let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(argument)
                    {
                        // Its complete index identity belongs to this port, including
                        // scoped labels. Count it in written order on the owning tensor.
                        result.add(slot, 1);
                        continue;
                    }
                    child.count_indices = count_arguments && Self::nested_index_argument(argument);
                    if child.validate || child.count_indices {
                        result.append(&Self::observe_indices_with_interfaces(
                            argument,
                            slots,
                            cache,
                            child,
                            exclude_tensor_leaf,
                            proofs,
                        )?);
                    }
                }
                if inputs + outputs > 0 && (inputs != 1 || outputs != 1) {
                    return Err(TensorInferenceError::invalid(format!(
                        "tensor factor `{symbol}` requires exactly one direct `in` and one direct `out` placeholder",
                    )));
                }
            }
            _ => {}
        }
        if let Some(representation) = scope.pair_representation {
            let selected = result
                .entries()
                .any(|(slot, count)| slot.rep().base().rep == representation.base() && count >= 2);
            result.selected_pair |= selected;
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bis;
    use crate::representations::Bispinor;
    use crate::tensor::SymbolicTensor;
    use crate::tensor::inference::tests::{INFERENCE_CALLS, INFERENCE_VISITS};
    use crate::test_support;
    use spenso::network::library::symbolic::ETS;
    use spenso::structure::dimension::Dimension;
    use spenso::structure::representation::{ExtendibleReps, Minkowski};
    use spenso::{bracket, mink, tensor_symbol, vector_symbol};
    use symbolica::atom::{Atom, FunctionBuilder};
    use symbolica::{function, symbol};

    #[test]
    fn initial_sum_producer_transports_its_boundary_to_fresh_generic_terms() {
        test_support::test_initialize();
        let mu = mink!(4, 99469);
        let left = function!(tensor_symbol!("initial_generic_left"), &mu);
        let right = function!(tensor_symbol!("initial_generic_right"), &mu);
        let declared = SymbolicTensor::infer(left.clone()).unwrap().structure;
        SCOPED_INTERFACE_INFERENCES.with(|count| count.set(0));
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let value = SymbolicTensor::from_produced_parts(left + right, declared).unwrap();
        assert_eq!(value.rank(), 1);
        assert_eq!(SCOPED_INTERFACE_INFERENCES.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn parser_observations_keep_competing_sum_orders_local() {
        test_support::test_initialize();
        let mu = mink!(4, 99731);
        let nu = mink!(4, 99733);
        let shared = function!(tensor_symbol!("parser_shared_term"), &mu, &nu);
        let first = &shared + function!(tensor_symbol!("parser_first_term"), &mu, &nu);
        let second = &shared + function!(tensor_symbol!("parser_second_term"), &mu, &nu);
        let forward = SymbolicTensor::infer(shared.clone()).unwrap().structure;
        let reverse =
            PartialStructure::from_logical_slots(forward.logical_slots().into_iter().rev());
        assert_ne!(forward, reverse);
        let source = TensorProofs::default();
        source.retain_boundary(first.as_view(), &forward);
        source.retain_boundary(second.as_view(), &reverse);

        for (sum, declaration) in [(&first, &forward), (&second, &reverse), (&first, &forward)] {
            let (history, private) =
                InterfaceInference::parser_occurrences(sum.as_view(), &source, None).unwrap();
            assert_eq!(history.entries().count(), 2);
            assert_eq!(
                private
                    .unwrap()
                    .retained_interface(shared.as_view())
                    .as_ref(),
                Some(declaration)
            );
            assert!(source.retained_interface(shared.as_view()).is_none());
        }
        // A term with its own declaration keeps that logical order in either sum.
        source.retain_boundary(shared.as_view(), &reverse);
        let (_, private) =
            InterfaceInference::parser_occurrences(first.as_view(), &source, None).unwrap();
        assert_eq!(
            private
                .unwrap()
                .retained_interface(shared.as_view())
                .as_ref(),
            Some(&reverse)
        );
        assert_eq!(source.retained_interface(shared.as_view()), Some(reverse));
    }

    #[test]
    fn provisional_parser_boundary_preserves_existing_raw_composite_ports() {
        test_support::test_initialize();
        let mu = mink!(4, 99735);
        let value = bracket!(
            function!(tensor_symbol!("parser_raw_composite_left"), &mu),
            function!(tensor_symbol!("parser_raw_composite_right"), &mu)
        );
        let source = TensorProofs::default();
        let raw = vec![mu.clone(), mu];
        let leaves = source.parser_leaf_interfaces();
        leaves
            .lock()
            .unwrap()
            .insert(value.clone(), (true, raw.clone()));
        let scalar = PartialStructure::from_logical_slots([]);
        let (_, private) =
            InterfaceInference::parser_occurrences(value.as_view(), &source, Some(&scalar))
                .unwrap();
        let private = private.unwrap();
        assert_eq!(private.retained_interface(value.as_view()), Some(scalar));
        assert_eq!(
            private.parser_leaf_interfaces().lock().unwrap()[&value].1,
            raw
        );
        assert_eq!(leaves.lock().unwrap()[&value].1, raw);
        assert!(source.retained_interface(value.as_view()).is_none());
    }

    #[test]
    fn initial_producer_scopes_compose_constructor_interfaces_without_inference() {
        test_support::test_initialize();
        let mu = mink!(4, 99465);
        let nu = mink!(4, 99467);
        let left = function!(vector_symbol!("initial_producer_p"), &mu)
            * function!(vector_symbol!("initial_producer_q"), &mu);
        let right = function!(vector_symbol!("initial_producer_r"), &nu)
            * function!(vector_symbol!("initial_producer_s"), &nu);
        SCOPED_INTERFACE_INFERENCES.with(|count| count.set(0));
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let value = SymbolicTensor::from_produced_parts(
            (left + right).pow(2),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        assert!(value.is_scalar());
        assert_eq!(SCOPED_INTERFACE_INFERENCES.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn producer_occurrence_scopes_reuse_admitted_interfaces() {
        test_support::test_initialize();
        let mu = mink!(4, 99461);
        let a = function!(tensor_symbol!("producer_occurrence_a"), &mu);
        let b = function!(tensor_symbol!("producer_occurrence_b"), &mu);
        let slot = SlotMatcher::default()
            .parse::<LibraryRep, AbstractIndex>(mu.as_view())
            .unwrap();
        for (expression, expected) in [(&a + &b, 1), (bracket!(&a, &b).pow(2), 0)] {
            let source = SymbolicTensor::infer(expression).unwrap();
            SCOPED_INTERFACE_INFERENCES.with(|count| count.set(0));
            INFERENCE_CALLS.with(|count| count.set(0));
            INFERENCE_VISITS.with(|count| count.set(0));
            let occurrences = InterfaceInference::produced_occurrences(
                source.expression().as_view(),
                &mut SlotMatcher::default(),
                &source.proofs,
            )
            .unwrap();
            assert_eq!(occurrences.count(slot), expected);
            assert_eq!(SCOPED_INTERFACE_INFERENCES.with(|count| count.get()), 0);
            assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
            assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        }
    }

    #[test]
    fn producer_occurrence_scopes_reject_missing_facts_without_inference() {
        test_support::test_initialize();
        let mu = mink!(4, 99463);
        let a = function!(tensor_symbol!("producer_missing_a"), &mu);
        let source = SymbolicTensor::infer(a.clone()).unwrap();
        let unknown = a + function!(tensor_symbol!("producer_missing_b"), &mu);
        SCOPED_INTERFACE_INFERENCES.with(|count| count.set(0));
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        assert!(
            InterfaceInference::produced_occurrences(
                unknown.as_view(),
                &mut SlotMatcher::default(),
                &source.proofs,
            )
            .is_err()
        );
        assert_eq!(SCOPED_INTERFACE_INFERENCES.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn one_observation_checks_scopes_and_keeps_branch_multiplicity() {
        test_support::test_initialize();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let a = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99401));
        let leaf = |name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(a.to_atom())
                .finish()
        };
        let alternatives = leaf("observation_A") + leaf("observation_B");
        let valid = &alternatives * leaf("observation_C");
        let invalid = valid.clone() * leaf("observation_D");
        for (value, expected) in [(alternatives, 1), (valid, 2), (invalid, 3)] {
            let occurrences = InterfaceInference::scoped_occurrences(
                value.as_view(),
                &mut SlotMatcher::default(),
            )
            .unwrap();
            assert_eq!(occurrences.count(a), expected);
            assert_eq!(occurrences.validate().is_ok(), expected <= 2);
        }
    }

    #[test]
    fn additive_branches_keep_internal_dummy_indices_local() {
        crate::test_support::test_initialize();
        let mu = spenso::mink!(4, 99451);
        let nu = spenso::mink!(4, 99453);
        let p = FunctionBuilder::new(spenso::vector_symbol!("add_scope_p"))
            .add_arg(&mu)
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("add_scope_q"))
            .add_arg(&mu)
            .finish();
        let r = FunctionBuilder::new(spenso::vector_symbol!("add_scope_r"))
            .add_arg(&nu)
            .finish();
        let metric = FunctionBuilder::new(ETS.metric)
            .add_args([&mu, &nu])
            .finish();
        let x = Atom::var(symbolica::symbol!("add_scope_x"));
        for expression in [
            (&x + p.pow(2)) * &p,
            (&x + &p * &q) * &p,
            (&metric * &p + &r) * &p,
        ] {
            SymbolicTensor::<PartialStructure>::validate_atom(&expression).unwrap();
            let value = SymbolicTensor::infer(expression.clone()).unwrap();
            assert_eq!(value.expression(), &expression);
            let occurrences = InterfaceInference::scoped_occurrences(
                expression.as_view(),
                &mut SlotMatcher::default(),
            )
            .unwrap();
            occurrences.validate().unwrap();
            assert!(
                InterfaceInference::has_explicit_index_pairs(
                    expression.as_view(),
                    Minkowski {}.into(),
                    |_| false,
                )
                .unwrap()
            );
        }
        // An explicitly supplied cubic term has a different index scope.
        // Construct it directly; do not expand a factored tensor expression.
        let cubic = &x * &p + p.pow(3);
        assert!(SymbolicTensor::infer(cubic).is_err());
        let invalid_branch = (&x + p.pow(3)) * &q;
        assert!(SymbolicTensor::<PartialStructure>::validate_atom(&invalid_branch).is_err());
        assert!(
            InterfaceInference::has_explicit_index_pairs(
                invalid_branch.as_view(),
                Minkowski {}.into(),
                |_| false,
            )
            .is_err()
        );
    }

    #[test]
    fn admitted_slot_labels_are_not_chain_placeholder_syntax() {
        test_support::test_initialize();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let head = tensor_symbol!("slot_scope_label_T");
        for label in [SPENSO_TAG.chain_in, SPENSO_TAG.chain_out] {
            for index in [
                AbstractIndex::Symbol(label.into()),
                AbstractIndex::Symbol(label.into()).scoped(symbol!("slot_scope_label")),
            ] {
                let slot = rep.slot(index);
                let expression = function!(head, slot.to_atom());
                let occurrences = InterfaceInference::scoped_occurrences(
                    expression.as_view(),
                    &mut SlotMatcher::default(),
                )
                .unwrap();
                assert_eq!(occurrences.entries().collect::<Vec<_>>(), [(slot, 1)]);
                let source = SymbolicTensor::infer(expression.clone()).unwrap();
                assert_eq!(
                    source.structure.logical_slots(),
                    vec![rep.slot(PartialIndex::Explicit(index))]
                );
                assert_eq!(
                    source.contract(Default::default()).unwrap().expression,
                    expression
                );
            }
            let compact = function!(head, Minkowski {}.new_rep(label).to_symbolic([]));
            SymbolicTensor::infer(compact).unwrap();
            let actual_placeholder = function!(head, label);
            assert!(
                InterfaceInference::scoped_occurrences(
                    actual_placeholder.as_view(),
                    &mut SlotMatcher::default(),
                )
                .is_err()
            );
        }
    }

    #[test]
    fn cached_scope_proof_is_local_and_scalar_metadata_is_opaque() {
        test_support::test_initialize();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let input = Atom::var(SPENSO_TAG.chain_in);
        let output = Atom::var(SPENSO_TAG.chain_out);
        let leaf = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("observation_channel"))
            .add_arg(&input)
            .add_arg(&output)
            .finish();
        let mut slots = SlotMatcher::default();
        let mut cache = Cache::new();
        let in_factor = ObservationScope {
            validate: true,
            inside_factor: true,
            count_indices: true,
            ..ObservationScope::default()
        };
        assert!(
            InterfaceInference::observe_indices(
                leaf.as_view(),
                &mut slots,
                &mut cache,
                in_factor,
                &mut |_| false
            )
            .is_ok()
        );
        assert!(
            InterfaceInference::observe_indices(
                leaf.as_view(),
                &mut slots,
                &mut cache,
                ObservationScope {
                    inside_factor: false,
                    ..in_factor
                },
                &mut |_| false,
            )
            .is_err()
        );
        let chain = FunctionBuilder::new(SPENSO_TAG.chain)
            .add_arg(rep.slot::<AbstractIndex, _>(99403).to_atom())
            .add_arg(rep.slot::<AbstractIndex, _>(99405).to_atom())
            .add_arg(leaf)
            .finish();
        assert!(InterfaceInference::scoped_occurrences(chain.as_view(), &mut slots).is_ok());
        let nested = FunctionBuilder::new(SPENSO_TAG.trace)
            .add_arg(rep.to_symbolic([]))
            .add_arg(chain)
            .finish();
        assert!(InterfaceInference::scoped_occurrences(nested.as_view(), &mut slots).is_err());
        let scalar = FunctionBuilder::new(symbol!("observation_metadata"; Scalar))
            .add_arg(nested)
            .add_arg(input)
            .finish();
        assert!(
            InterfaceInference::scoped_occurrences(scalar.as_view(), &mut slots)
                .unwrap()
                .is_empty()
        );
    }
    #[test]
    fn scalar_inverse_scopes_keep_invalid_local_indices_and_placeholders_visible() {
        test_support::test_initialize();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let slot = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99417));
        let leaf = |name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(slot.to_atom())
                .finish()
        };
        let a = leaf("inverse_scope_A");
        let b = leaf("inverse_scope_B");
        let c = leaf("inverse_scope_C");
        for denominator in [&a * &b, a.pow(2), &a * &b + Atom::one()] {
            let inverse = bracket!(denominator).pow(-1);
            let value = &a * inverse;
            let occurrences = InterfaceInference::scoped_occurrences(
                value.as_view(),
                &mut SlotMatcher::default(),
            )
            .unwrap();
            assert_eq!(occurrences.count(slot), 1);
            occurrences.validate().unwrap();
        }
        let invalid = bracket!(&a * &b * &c).pow(-1);
        let occurrences =
            InterfaceInference::scoped_occurrences(invalid.as_view(), &mut SlotMatcher::default())
                .unwrap();
        assert!(occurrences.validate().is_err());
        let malformed = bracket!(&a * &b * Atom::var(SPENSO_TAG.chain_in)).pow(-1);
        assert!(
            InterfaceInference::scoped_occurrences(
                malformed.as_view(),
                &mut SlotMatcher::default()
            )
            .is_err()
        );
    }

    #[test]
    fn powers_share_closed_dummy_scopes_and_open_port_multiplicity() {
        test_support::test_initialize();
        let slot = mink!(4, 99431);
        let leaf = |name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(&slot)
                .finish()
        };
        let a = leaf("power_scope_A");
        let b = leaf("power_scope_B");
        let outer = leaf("power_scope_outer");
        let closed = bracket!(&a * &b);
        for exponent in [-2, -1, 2, 7] {
            let expression = closed.pow(exponent) * &outer;
            SymbolicTensor::<PartialStructure>::validate_atom(&expression).unwrap();
            let value = SymbolicTensor::infer(expression.clone()).unwrap();
            assert_eq!(value.expression(), &expression);
            assert_eq!(value.structure().canonical().external_structure().len(), 1);
            let counts = InterfaceInference::index_occurrences(
                expression.as_view(),
                &mut SlotMatcher::default(),
                &mut Cache::default(),
            );
            assert_eq!(
                counts.entries().map(|(_, count)| count).collect::<Vec<_>>(),
                vec![1]
            );
            assert!(
                InterfaceInference::has_explicit_index_pairs(
                    expression.as_view(),
                    Minkowski {}.into(),
                    |_| false
                )
                .unwrap()
            );
        }
        let invalid = a.pow(2) * &outer;
        assert!(SymbolicTensor::<PartialStructure>::validate_atom(&invalid).is_err());
        assert!(SymbolicTensor::infer(invalid.clone()).is_err());
        assert!(
            InterfaceInference::has_explicit_index_pairs(
                invalid.as_view(),
                Minkowski {}.into(),
                |_| false
            )
            .is_err()
        );
        assert!(
            InterfaceInference::index_occurrences(
                invalid.as_view(),
                &mut SlotMatcher::default(),
                &mut Cache::default()
            )
            .validate()
            .is_err()
        );
    }

    #[test]
    fn closed_exponents_keep_their_contractions_independent_of_outer_indices() {
        test_support::test_initialize();
        let mu = mink!(4, 99621);
        let nu = mink!(4, 99623);
        let leaf = |name, slot: &Atom| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(slot)
                .finish()
        };
        let c = leaf("exponent_scope_C", &mu);
        let d = leaf("exponent_scope_D", &mu);
        for inner in [&mu, &nu] {
            let exponent = leaf("exponent_scope_A", inner) * leaf("exponent_scope_B", inner);
            let power = Atom::var(Symbol::E).pow(&exponent);
            for (outer, count, rank) in [(c.clone(), 1, 1), (&c * &d, 2, 0)] {
                let expression = &power * outer;
                let occurrences = InterfaceInference::scoped_occurrences(
                    expression.as_view(),
                    &mut SlotMatcher::default(),
                )
                .unwrap();
                assert_eq!(
                    occurrences
                        .entries()
                        .map(|(_, count)| count)
                        .collect::<Vec<_>>(),
                    vec![count]
                );
                occurrences.validate().unwrap();
                let value = SymbolicTensor::infer(expression.clone()).unwrap();
                assert_eq!(value.expression(), &expression);
                assert_eq!(
                    value.structure().canonical().external_structure().len(),
                    rank
                );
                SymbolicTensor::<PartialStructure>::validate_atom(&expression).unwrap();
            }
            assert!(
                InterfaceInference::scoped_occurrences(
                    power.as_view(),
                    &mut SlotMatcher::default(),
                )
                .unwrap()
                .is_empty()
            );
            assert!(
                InterfaceInference::has_explicit_index_pairs(
                    power.as_view(),
                    Minkowski {}.into(),
                    |_| false,
                )
                .unwrap()
            );
        }
        let metadata = symbol!("exponent_scope_metadata"; Scalar);
        let exponent = FunctionBuilder::new(metadata).add_arg(&c * &d).finish();
        let power = Atom::var(Symbol::E).pow(exponent);
        let expression = &power * &c * &d;
        assert_eq!(
            SymbolicTensor::infer(expression.clone())
                .unwrap()
                .expression(),
            &expression
        );
        assert!(
            !InterfaceInference::has_explicit_index_pairs(
                power.as_view(),
                Minkowski {}.into(),
                |_| false,
            )
            .unwrap()
        );
    }

    #[test]
    fn nested_scalar_exponents_retain_pair_queries_without_exporting_indices() {
        test_support::test_initialize();
        let mu = mink!(4, 99625);
        let leaf = |name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(&mu)
                .finish()
        };
        let closed = leaf("nested_exponent_A") * leaf("nested_exponent_B");
        let metadata = FunctionBuilder::new(symbol!("nested_exponent_metadata"; Scalar))
            .add_arg(&closed)
            .finish();
        let base = Atom::var(symbol!("nested_exponent_base"));
        for (mut exponent, expected_pairs) in [(closed, true), (metadata, false)] {
            for depth in 1..=3 {
                exponent = base.pow(exponent);
                assert!(
                    InterfaceInference::scoped_occurrences(
                        exponent.as_view(),
                        &mut SlotMatcher::default(),
                    )
                    .unwrap()
                    .is_empty(),
                    "depth {depth} must not export closed exponent indices",
                );
                let value = SymbolicTensor::infer(exponent.clone()).unwrap();
                assert_eq!(value.expression(), &exponent);
                assert!(value.structure().canonical().is_scalar());
                assert_eq!(
                    InterfaceInference::has_explicit_index_pairs(
                        exponent.as_view(),
                        Minkowski {}.into(),
                        |_| false,
                    )
                    .unwrap(),
                    expected_pairs,
                    "depth {depth} must retain nested pairs while respecting Scalar metadata",
                );
                assert!(
                    !InterfaceInference::has_explicit_index_pairs(
                        exponent.as_view(),
                        Minkowski {}.into(),
                        |_| true,
                    )
                    .unwrap(),
                    "depth {depth} must respect excluded tensor leaves",
                );
                assert!(
                    !InterfaceInference::has_explicit_index_pairs(
                        exponent.as_view(),
                        Bispinor {}.into(),
                        |_| false,
                    )
                    .unwrap(),
                    "depth {depth} must match the requested representation",
                );
            }
        }
    }

    #[test]
    fn exponent_scoping_cannot_hide_open_ports_or_invalid_local_multiplicity() {
        test_support::test_initialize();
        let mu = mink!(4, 99631);
        let leaf = |name, slot| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(slot)
                .finish()
        };
        let a = leaf("exponent_invalid_A", mu.clone());
        let b = leaf("exponent_invalid_B", mu.clone());
        let c = leaf("exponent_invalid_C", mu);
        let open = leaf("exponent_invalid_open", mink!(4));
        for exponent in [
            a.clone(),
            open,
            &a * &b * &c,
            Atom::var(SPENSO_TAG.chain_in),
        ] {
            let expression = Atom::var(Symbol::E).pow(exponent);
            assert!(SymbolicTensor::<PartialStructure>::validate_atom(&expression).is_err());
            if InterfaceInference::contains_tensor_symbol(expression.as_view()) {
                assert!(SymbolicTensor::infer(expression.clone()).is_err());
            } else {
                assert!(
                    SymbolicTensor::infer(expression.clone())
                        .unwrap()
                        .is_scalar()
                );
            }
            assert!(
                InterfaceInference::has_explicit_index_pairs(
                    expression.as_view(),
                    Minkowski {}.into(),
                    |_| true,
                )
                .is_err()
            );
        }
        let invalid = Atom::var(Symbol::E).pow(&a * &b * &c);
        assert!(
            InterfaceInference::index_occurrences(
                invalid.as_view(),
                &mut SlotMatcher::default(),
                &mut Cache::default(),
            )
            .validate()
            .is_err()
        );
    }

    #[test]
    fn exponents_have_independent_chain_placeholder_scopes() {
        test_support::test_initialize();
        let matrix = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("exponent_scope_matrix"))
            .add_args([SPENSO_TAG.chain_in, SPENSO_TAG.chain_out])
            .finish();
        let trace = FunctionBuilder::new(SPENSO_TAG.trace)
            .add_arg(Bispinor {}.new_rep(4).to_symbolic([]))
            .add_arg(&matrix)
            .finish();
        let chain = |exponent| {
            FunctionBuilder::new(SPENSO_TAG.chain)
                .add_args([bis!(4, 99641), bis!(4, 99643)])
                .add_arg(Atom::var(Symbol::E).pow(exponent) * &matrix)
                .finish()
        };
        let valid = chain(trace);
        let value = SymbolicTensor::infer(valid.clone()).unwrap();
        assert_eq!(value.structure().canonical().external_structure().len(), 2);
        SymbolicTensor::<PartialStructure>::validate_atom(&valid).unwrap();
        // The same tensor ports are legal as outer matrix factors, but cannot
        // borrow that scope when used as a naked exponent.
        let malformed = FunctionBuilder::new(SPENSO_TAG.trace)
            .add_arg(Bispinor {}.new_rep(4).to_symbolic([]))
            .add_arg(
                FunctionBuilder::new(SPENSO_TAG.tensor_symbol("exponent_scope_bad_matrix"))
                    .add_args([SPENSO_TAG.chain_in, SPENSO_TAG.chain_in])
                    .finish(),
            )
            .finish();
        for exponent in [matrix.clone(), malformed] {
            let invalid = chain(exponent);
            assert!(SymbolicTensor::infer(invalid.clone()).is_err());
            assert!(SymbolicTensor::<PartialStructure>::validate_atom(&invalid).is_err());
        }
    }

    #[test]
    fn powered_composites_repeat_only_boundary_ports() {
        crate::test_support::test_initialize();
        let mu = spenso::mink!(4, 99441);
        let nu = spenso::mink!(4, 99443);
        let p = FunctionBuilder::new(spenso::vector_symbol!("power_boundary_p"))
            .add_arg(&mu)
            .finish();
        let outer = FunctionBuilder::new(spenso::vector_symbol!("power_boundary_outer"))
            .add_arg(&mu)
            .finish();
        let metric = FunctionBuilder::new(ETS.metric)
            .add_arg(&mu)
            .add_arg(&nu)
            .finish();
        let base = spenso::bracket!(metric * &p);
        let x = Atom::var(symbolica::symbol!("power_boundary_x"));
        let unbracketed = FunctionBuilder::new(ETS.metric)
            .add_arg(&mu)
            .add_arg(&nu)
            .finish()
            * &p;
        for (base, normalized_base) in [
            (base.clone(), base.clone()),
            (&base * &x + &base, &unbracketed * &x + &unbracketed),
        ] {
            let expression = base.pow(2) * &outer;
            SymbolicTensor::<PartialStructure>::validate_atom(&expression).unwrap();
            let value = SymbolicTensor::infer(expression.clone()).unwrap();
            assert_eq!(value.expression(), &(normalized_base.pow(2) * &outer));
            assert_eq!(value.structure().canonical().external_structure().len(), 1);
            assert!(
                InterfaceInference::has_explicit_index_pairs(
                    expression.as_view(),
                    Minkowski {}.into(),
                    |_| false
                )
                .unwrap()
            );
            let contracted = value
                .contract(Default::default())
                .unwrap()
                .to_dots()
                .unwrap();
            SymbolicTensor::<PartialStructure>::validate_atom(contracted.expression()).unwrap();
            assert_eq!(contracted.structure(), value.structure());
        }
    }

    #[test]
    fn explicit_pair_query_preserves_factored_branches_and_exclusions() {
        crate::test_support::test_initialize();
        let mu = spenso::mink!(4, 99501);
        let nu = spenso::mink!(4, 99503);
        let leaf = |name, slot: &Atom| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(slot)
                .finish()
        };
        let a = leaf("incidence_A", &mu);
        let b = leaf("incidence_B", &mu);
        let c = leaf("incidence_C", &mu);
        let d = leaf("incidence_D", &nu);
        let x = Atom::var(symbolica::symbol!("incidence_x"));
        let y = Atom::var(symbolica::symbol!("incidence_y"));
        let pairs = |value: &Atom| {
            InterfaceInference::has_explicit_index_pairs(
                value.as_view(),
                Minkowski {}.into(),
                |_| false,
            )
        };
        assert!(!pairs(&(&a + &b)).unwrap());
        assert!(pairs(&((&a + &b) * &c)).unwrap());
        assert!(!pairs(&((&a + Atom::one()) * &d)).unwrap());
        assert!(pairs(&((&a + Atom::one()) * &b)).unwrap());
        assert!(pairs(&((&a + &b) * &b * &c)).is_err());
        let weighted = (&x - &y) * &a * &b * (&x + &y).pow(30);
        let original = weighted.clone();
        assert!(pairs(&weighted).unwrap());
        assert_eq!(weighted, original);

        let momentum = SPENSO_TAG.tensor_symbol("incidence_momentum");
        let p = FunctionBuilder::new(momentum).add_arg(&mu).finish();
        assert!(
            !InterfaceInference::has_explicit_index_pairs(
                (&a * &p).as_view(),
                Minkowski {}.into(),
                |value| matches!(value, AtomView::Fun(function) if function.get_symbol() == momentum),
            )
            .unwrap()
        );
        // Exclusions cannot hide a malformed original three-port scope.
        let q = FunctionBuilder::new(momentum)
            .add_arg(1)
            .add_arg(&mu)
            .finish();
        let malformed = (&p - &q) * &a * &b;
        assert!(
            InterfaceInference::has_explicit_index_pairs(
                malformed.as_view(),
                Minkowski {}.into(),
                |value| matches!(value,AtomView::Fun(function) if function.get_symbol()==momentum),
            )
            .unwrap_err()
            .to_string()
            .contains("more than two compatible tensor ports")
        );
        // Equalized momentum leaves must not cancel two distinct branches.
        // Their independent nu port keeps this witness well formed before
        // excluding momenta; A(mu)*B(mu) is the pair being detected.
        let free_p = FunctionBuilder::new(momentum).add_arg(&nu).finish();
        let free_q = FunctionBuilder::new(momentum)
            .add_arg(1)
            .add_arg(&nu)
            .finish();
        let weighted = (&free_p - &free_q) * &a * &b;
        assert!(
            InterfaceInference::has_explicit_index_pairs(
                weighted.as_view(),
                Minkowski {}.into(),
                |value| matches!(value,AtomView::Fun(function) if function.get_symbol()==momentum),
            )
            .unwrap()
        );
        // The diagnostic query also accepts nested alternatives whose
        // interfaces differ. Failed scope inference must leave their original
        // occurrences available, rather than introduce typed admission here.
        let nested = ((&a + Atom::one()) * &d + &d) * &b;
        assert!(pairs(&nested).unwrap());
        let nested_filtered = ((&a * &p + Atom::one()) * &d + &d) * &b;
        assert!(
            !InterfaceInference::has_explicit_index_pairs(
                nested_filtered.as_view(),
                Minkowski {}.into(),
                |value| matches!(value, AtomView::Fun(function) if function.get_symbol() == momentum),
            )
            .unwrap()
        );
        let scalar = FunctionBuilder::new(symbolica::symbol!("incidence_metadata"; Scalar))
            .add_arg(&b * &c)
            .finish();
        assert!(!pairs(&(&a * scalar)).unwrap());
        assert!(
            InterfaceInference::has_explicit_index_pairs(
                spenso::bracket!(&a * &b).as_view(),
                Minkowski {}.into(),
                |value| matches!(value, AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.bracket),
            )
            .unwrap()
        );
    }

    #[test]
    fn explicit_pair_query_keeps_power_scopes_and_rejects_malformed_indices() {
        crate::test_support::test_initialize();
        let mu = spenso::mink!(4, 99511);
        let leaf = |name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(&mu)
                .finish()
        };
        let a = leaf("incidence_power_A");
        let b = leaf("incidence_power_B");
        let c = leaf("incidence_power_C");
        let pairs = |value: &Atom| {
            InterfaceInference::has_explicit_index_pairs(
                value.as_view(),
                Minkowski {}.into(),
                |_| false,
            )
        };
        assert!(pairs(&a.pow(2)).unwrap());
        assert!(pairs(&a.pow(3)).is_err());
        let scalar = spenso::bracket!(&a * &b);
        for power in [-2, -1, 2, 7] {
            assert!(pairs(&(scalar.pow(power) * &c)).unwrap());
        }
        assert!(pairs(&spenso::bracket!(&a * &b * &c).pow(-1)).is_err());
        assert!(pairs(&Atom::var(SPENSO_TAG.chain_in)).is_err());
        let momentum = SPENSO_TAG.tensor_symbol("incidence_power_momentum");
        let p = FunctionBuilder::new(momentum).add_arg(&mu).finish();
        for power in [-1, 2, 7] {
            let value = spenso::bracket!(&a * &p).pow(power) * &b;
            assert!(
                !InterfaceInference::has_explicit_index_pairs(
                    value.as_view(),
                    Minkowski {}.into(),
                    |value| matches!(value, AtomView::Fun(function) if function.get_symbol() == momentum),
                )
                .unwrap()
            );
        }
    }

    #[test]
    fn pair_exclusion_predicate_distinguishes_leaves_of_the_same_head() {
        let mu = spenso::mink!(4, 99581);
        let a = FunctionBuilder::new(spenso::tensor_symbol!("pair_leaf_filter::A"))
            .add_arg(&mu)
            .finish();
        let head = spenso::tensor_symbol!("pair_leaf_filter::v");
        let selected = FunctionBuilder::new(head).add_arg(0).add_arg(&mu).finish();
        let retained = FunctionBuilder::new(head).add_arg(1).add_arg(&mu).finish();
        let exclude = |value: AtomView<'_>| {
            matches!(value,AtomView::Fun(function)
            if function.get_symbol()==head && function.get(0).is_zero())
        };
        for (input, expected) in [
            (&a * &selected, false),
            (&a * &retained, true),
            (&a * (&selected + &retained), true),
            (spenso::bracket!(&a * &selected).pow(3), false),
            (spenso::bracket!(&a * &retained).pow(-2), true),
        ] {
            let original = input.clone();
            assert_eq!(
                InterfaceInference::has_explicit_index_pairs(
                    input.as_view(),
                    Minkowski {}.into(),
                    exclude
                )
                .unwrap(),
                expected
            );
            assert_eq!(input, original);
        }
    }

    #[test]
    fn pair_exclusion_predicate_does_not_replay_leaf_normalizers() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static CALLS: AtomicUsize = AtomicUsize::new(0);
        let mu = spenso::mink!(4, 99591);
        let head = spenso::tensor_symbol!(
            "pair_leaf_filter::callback",
            norm = |_, _| {
                CALLS.fetch_add(1, Ordering::SeqCst);
            }
        );
        let left = FunctionBuilder::new(head).add_arg(0).add_arg(&mu).finish();
        let right = FunctionBuilder::new(head).add_arg(1).add_arg(&mu).finish();
        let input = &left * &right;
        let original = input.clone();
        let calls = CALLS.load(Ordering::SeqCst);
        assert!(!InterfaceInference::has_explicit_index_pairs(input.as_view(),Minkowski {}.into(),
            |value| matches!(value,AtomView::Fun(function) if function.get_symbol()==head && function.get(0).is_zero())
        ).unwrap());
        assert_eq!(CALLS.load(Ordering::SeqCst), calls);
        assert_eq!(input, original);
    }
    #[test]
    fn pair_exclusions_cannot_hide_original_three_port_multiplicity() {
        let mu = spenso::mink!(4, 99597);
        let a = FunctionBuilder::new(spenso::tensor_symbol!("pair_scope_validation::A"))
            .add_arg(&mu)
            .finish();
        let head = spenso::tensor_symbol!("pair_scope_validation::excluded");
        let first = FunctionBuilder::new(head).add_arg(1).add_arg(&mu).finish();
        let second = FunctionBuilder::new(head).add_arg(2).add_arg(&mu).finish();
        let source = &a * &first * &second;
        let original = source.clone();
        let error = InterfaceInference::has_explicit_index_pairs(
            source.as_view(),
            Minkowski {}.into(),
            |value| matches!(value,AtomView::Fun(function) if function.get_symbol()==head),
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("more than two compatible tensor ports")
        );
        assert_eq!(source, original);
        assert!(
            !InterfaceInference::has_explicit_index_pairs(
                (&a * &first).as_view(),
                Minkowski {}.into(),
                |value| matches!(value,AtomView::Fun(function) if function.get_symbol()==head)
            )
            .unwrap()
        );
    }
}
