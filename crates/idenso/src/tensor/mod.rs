use self::composition::ExplicitIndexOccurrences;
use self::inference::{InterfaceInference, LeafInference, TensorInferenceError, TensorInterface};
use self::raw_leaves::RawLeafInterfaces;
use self::regions::AdmittedInterfaces;
use crate::NetworkToolingError;
use delegate::delegate;
use eyre::eyre;
use linnet::half_edge::subgraph::subset::SubSet;
use spenso::algebra::ScalarMul;
use spenso::bracket;
use spenso::contraction::{Contract, ContractionError, Trace};
use spenso::iterators::IteratableTensor;
use spenso::network::library::function_lib::Wrap;
use spenso::network::library::symbolic::{ETS, ExplicitKey, TensorLibrary};
use spenso::network::library::{DummyKey, DummyLibrary, FunctionLibrary, FunctionLibraryError};
use spenso::network::parsing::{
    AtomStructureExt, LeafInterfaceCache, ParseSettings, ShadowedStructure, StrictTensorFilter,
    StructureFromAtom, TensorFromExpression, TensorLibraryFor,
};
use spenso::network::store::NetworkStore;
use spenso::network::tags::SPENSO_TAG;
use spenso::network::{
    ExecutionResult, Network, Ref, Sequential, SmallestDegree, TensorNetworkError,
};
use spenso::shadowing::Concretize;
use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use spenso::structure::abstract_index::{AIND_SYMBOLS, AbstractIndex};
use spenso::structure::concrete_index::{ExpandedIndex, FlatIndex};
use spenso::structure::dimension::Dimension;
use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::representation::{LibraryRep, LibrarySlot, Representation};
use spenso::structure::slot::SlotError;
use spenso::structure::slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind, SlotMatcher};
use spenso::structure::{
    ApplyPendingIndexPermutation, CanonicalLayout, Canonicalized, HasName, HasStructure,
    IndexedStructure, MergeInfo, NamedStructure, OrderedStructure, PendingIndexPermutation,
    Reindexed, ScalarStructure, ScalarTensor, SlotIndex, StructureContract, StructureError,
    TensorIdentity, TensorShell, TensorStructure, ToSymbolic,
};
use spenso::tensors::parametric::MixedTensor;
use std::collections::{HashMap, HashSet};
use std::ops::AddAssign;
use std::sync::{Arc, OnceLock};
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol};
use symbolica::function;
use symbolica_utils::{AtomPrintExt, IntoArgs, IntoSymbol};

mod canonicalize;
mod collection;
pub(crate) mod planning;
#[cfg(test)]
pub(crate) use collection::SELECTED_OPENS;
pub(crate) use collection::{CollectedOutput, CollectedPolynomial, CollectionMode};
pub mod composition;
mod indexed;
pub mod inference;
mod network;
#[cfg(feature = "bincode")]
mod persistence;
mod raw_leaves;
mod regions;
mod replacement;
pub mod simplification;
pub use replacement::TensorRule;

#[cfg(test)]
mod network_metadata_tests;
#[cfg(test)]
pub mod tests;

impl<S: TensorInterface> FunctionLibrary<SymbolicTensor<S>, Atom> for Wrap {
    type Key = Symbol;
    fn apply(
        &self,
        key: &Self::Key,
        tensor: SymbolicTensor<S>,
    ) -> eyre::Result<SymbolicTensor<S>, FunctionLibraryError<Self::Key>> {
        let expression = function!(*key, &tensor.expression);
        #[cfg(debug_assertions)]
        if !InterfaceInference::intrinsic_normalization_head(*key)
            || !InterfaceInference::normalization_is_intrinsic(tensor.expression.as_view())
        {
            let structure = tensor.structure.to_partial_interface(&expression);
            TensorProofs::assert_normalization_preserves_interface(
                "Wrap::apply",
                &[&tensor.expression],
                &expression,
                structure.as_ref(),
                Some(*key),
            );
        }
        let interface = tensor
            .proofs
            .retained_interface(tensor.expression.as_view());
        let mut proofs = tensor.proofs.clone();
        if let Some(interface) = interface {
            proofs.publish_interfaces(expression.as_view(), interface, []);
        }
        Ok(SymbolicTensor {
            work: Default::default(),
            proofs,
            structure: tensor.structure.into_algebra_result(),
            is_composite: true,
            is_metric: false,
            expression,
        })
    }

    fn apply_scalar(
        &self,
        key: &Self::Key,
        scalar: Atom,
    ) -> eyre::Result<Atom, FunctionLibraryError<Self::Key>> {
        Ok(function!(*key, scalar))
    }
}

impl<Aind: ParseableAind + AbsInd + DummyAind> StructureFromAtom
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
    fn structure_from_atom(
        value: AtomView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        let structure = OrderedStructure::<LibraryRep, Aind>::structure_from_atom(value, matcher)?;
        let proofs = TensorProofs::parsed_leaf(value, &structure, None)
            .map_err(|error| StructureError::ParsingError(error.to_string()))?;
        Ok(structure.map_canonical(|structure| {
            let (is_composite, is_metric) = if let AtomView::Fun(f) = value {
                (false, f.get_symbol() == ETS.metric)
            } else {
                (true, false)
            };
            SymbolicTensor {
                work: Default::default(),
                proofs,
                structure,
                is_composite,
                is_metric,
                expression: value.to_owned(),
            }
        }))
    }
}

impl<'src, Aind, E, Sc, K, Lib, FunLib>
    TensorFromExpression<'src, OrderedStructure<LibraryRep, Aind>, Sc, K, Symbol, Aind, Lib, FunLib>
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>, E>
where
    Aind: AbsInd + DummyAind + ParseableAind,
    E: AtomCore + From<E::Output> + From<Atom> + From<AtomView<'src>>,
    K: std::fmt::Display,
    Lib: TensorLibraryFor<OrderedStructure<LibraryRep, Aind>, Self, Key = K>,
    FunLib: FunctionLibrary<Self, Sc, Key = Symbol>,
{
    type Admission = TensorProofs;

    fn admitted_structure(
        expression: AtomView<'_>,
        admission: &Self::Admission,
        next_index: impl FnMut() -> Aind,
    ) -> Result<Canonicalized<OrderedStructure<LibraryRep, Aind>>, StructureError> {
        admission.parser_structure(expression, &mut SlotMatcher::default(), next_index)
    }

    fn admitted_structure_with_matcher(
        expression: AtomView<'_>,
        admission: &Self::Admission,
        matcher: &mut SlotMatcher,
        next_index: impl FnMut() -> Aind,
    ) -> Result<Canonicalized<OrderedStructure<LibraryRep, Aind>>, StructureError> {
        admission.parser_structure(expression, matcher, next_index)
    }

    fn tensor_from_expression(
        expression: AtomOrView<'src>,
        structure: Canonicalized<OrderedStructure<LibraryRep, Aind>>,
        _tensor_library: &Lib,
        _function_library: &FunLib,
        settings: &ParseSettings,
        admission: Option<&Self::Admission>,
        _leaf_interfaces: Option<&LeafInterfaceCache>,
    ) -> Result<Canonicalized<Self>, TensorNetworkError<K, Symbol>> {
        let raw;
        let admission = if admission.is_none()
            && !matches!(expression.as_view(), AtomView::Fun(_))
            && !TensorProofs::uses_custom_graph_indices(&structure)
        {
            raw = TensorProofs::raw_parser_regions(
                expression.as_view(),
                settings.strict_tensor_filter,
            )
            .map_err(|error| TensorNetworkError::Other(eyre!(error)))?;
            Some(&raw)
        } else {
            admission
        };
        <Self as TensorFromExpression<
            'src,
            OrderedStructure<LibraryRep, Aind>,
            Sc,
            K,
            Symbol,
            Aind,
            Lib,
            FunLib,
        >>::tensor_from_leaf(expression, structure, admission)
    }

    fn tensor_from_leaf(
        expression: AtomOrView<'src>,
        structure: Canonicalized<OrderedStructure<LibraryRep, Aind>>,
        admission: Option<&Self::Admission>,
    ) -> Result<Canonicalized<Self>, TensorNetworkError<K, Symbol>> {
        let proofs = TensorProofs::parsed_leaf(expression.as_view(), &structure, admission)
            .map_err(|error| TensorNetworkError::Other(eyre!(error)))?;
        let (is_composite, is_metric) = if let AtomView::Fun(fun) = expression.as_view() {
            (false, fun.get_symbol() == ETS.metric)
        } else {
            (true, false)
        };
        let expression = match expression {
            AtomOrView::View(view) => view.into(),
            owned => owned.into_owned().into(),
        };
        Ok(structure.map_canonical(|structure| Self {
            work: Default::default(),
            proofs,
            structure,
            expression,
            is_metric,
            is_composite,
        }))
    }
}

/// A symbolic expression together with its tensor structure.
///
/// The default ordered, explicit-index structure supports network execution.
/// A [`spenso::structure::partial::PartialStructure`] also retains logical port
/// order and occurrence-local unresolved indices for positional composition.
/// The structure is independent of the expression, so a zero can retain its rank.
/// Tensor metadata is authoritative after admission. Algebra derives the
/// result's ports and scoped index uses from its operands; simplification of
/// the symbolic payload does not erase that construction history. Raw inputs
/// and arbitrary replacements must use a checked admission boundary. A callback
/// used by ordinary algebra must preserve the admitted tensor meaning.
/// The payload defaults to `Atom`; borrowed views use the same tensor carrier.
///
/// # Example
///
/// Admission separates a tensor's public ports from its scalar coefficients.
/// In particular, cancellation retains the rank and later contractions remain typed.
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
///
/// let port = Euclidean {}
///     .new_rep(3)
///     .to_lib()
///     .slot::<AbstractIndex, _>(AbstractIndex::Normal(7));
/// let vector = SymbolicTensor::infer(function!(tensor_symbol!("docs::vector"), port.to_atom()))?;
/// let zero = vector.try_sub(&vector)?;
/// assert!(zero.is_zero());
/// assert_eq!(zero.rank(), 1);
/// assert!(zero.multiply(&vector)?.is_scalar());
/// assert_eq!(vector.try_add(&vector)?.rank(), 1);
/// let two = SymbolicTensor::infer(Atom::num(2))?;
/// assert!(vector.pow(&two)?.is_scalar());
/// assert_eq!(vector.multiply(&two)?.rank(), 1);
///
/// // Raw parts must be admitted before they are used as authoritative metadata.
/// let admitted = SymbolicTensor::from_declared_parts(
///     vector.expression().clone(), vector.structure().clone(),
/// )?;
/// assert!(admitted.same_value(&vector));
/// let nominal_zero = SymbolicTensor::from_produced_parts(Atom::Zero, vector.structure().clone())?;
/// assert_eq!(nominal_zero.finish_parts()?.rank(), 1);
/// # Ok::<(), TensorInferenceError>(())
/// ```
///
/// Value comparisons exclude provenance and retained index-occurrence history.
/// `same_value` and `Hash` compare the symbolic value and logical port order;
/// `PartialEq` also compares the stored structure and classification flags.
/// Equal values can therefore carry different collision histories or completion
/// proofs. Consumers reusing such facts must compare the established facts too,
/// rather than treating value equality as a metadata certificate.
///
/// With the `bincode` feature, tensor archives are trusted GammaLoop state files.
/// Decode restores serialized interfaces, regional facts, and index history
/// without inference or replaying normalization in the same Symbolica state.
/// Symbolica's Atom decoder is not safe for untrusted bytes; the
/// archive magic and version distinguish supported state files, not a security boundary.
#[derive(Debug, Clone)]
pub struct SymbolicTensor<S = OrderedStructure<LibraryRep, AbstractIndex>, E = Atom> {
    pub(crate) structure: S,
    pub(crate) is_metric: bool,
    pub(crate) is_composite: bool,
    pub(crate) expression: E,
    pub(crate) proofs: TensorProofs,
    pub(crate) work: simplification::ReductionCache,
}

/// Admitted indices, regional interfaces, and the contracts for rewriting this payload.
/// This opaque type supplies parser admission context.
#[derive(Clone, Default)]
pub struct TensorProofs {
    /// Tensor product finishing is inert and there are no unresolved factor
    /// occurrences which ordinary multiplication could reorder or coalesce.
    pub(super) normalized_product: OnceLock<bool>,
    pub(crate) occurrences: OnceLock<Arc<ExplicitIndexOccurrences>>,
    /// Exact admitted subexpressions keep their free interfaces across known
    /// transformations; changed expressions cannot reuse another key's facts.
    pub(crate) admitted_interfaces: OnceLock<Arc<AdmittedInterfaces>>,
    pub(super) leaf_interfaces: RawLeafInterfaces,
}

impl std::fmt::Debug for TensorProofs {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("TensorProofs")
            .field("admitted", &self.is_admitted())
            .field("normalized_product", &self.normalized_product)
            .finish()
    }
}

/// The default symbolic network already stores AbstractIndex directly. Keep
/// its payload opaque when crossing between graph and semantic structures.
fn convert_parser_index<From: ParseableAind + 'static, To: ParseableAind + Copy + 'static>(
    index: From,
) -> Result<To, SlotError> {
    if let Some(&identity) = (&index as &dyn std::any::Any).downcast_ref::<To>() {
        Ok(identity)
    } else {
        To::from_view(index.to_atom().as_view()).map_err(Into::into)
    }
}

impl TensorProofs {
    /// Publication supplies index history even for a scalar or an empty tensor.
    pub(crate) fn is_admitted(&self) -> bool {
        self.occurrences.get().is_some()
    }

    fn parser_structure<Aind: AbsInd + ParseableAind>(
        &self,
        expression: AtomView<'_>,
        matcher: &mut SlotMatcher,
        mut next_index: impl FnMut() -> Aind,
    ) -> Result<Canonicalized<OrderedStructure<LibraryRep, Aind>>, StructureError> {
        let compact_metric = matches!(expression, AtomView::Fun(function)
            if InterfaceInference::is_compound_metric(function, matcher));
        // A compact metric is a partial contraction. Its free spectator axes
        // come from the regional interface, even if it has no direct raw slots.
        if !compact_metric
            && let Some(slots) = self
                .leaf_interfaces
                .with_fact(expression, |tensorial, ports| {
                    tensorial.then(|| {
                        ports
                            .map(|port| {
                                matcher.parse::<LibraryRep, Aind>(port).or_else(|_| {
                                    matcher
                                        .parse_representation::<LibraryRep>(port)
                                        .map(|representation| representation.slot(next_index()))
                                })
                            })
                            .collect::<Result<Vec<_>, _>>()
                    })
                })
                .flatten()
        {
            return Ok(OrderedStructure::new(slots?));
        }
        if let AtomView::Fun(function) = expression
            && !composition::is_composite_head(function.get_symbol(), &SPENSO_TAG)
            && function.get_symbol() != SPENSO_TAG.dot
            && !compact_metric
        {
            return Err(StructureError::ParsingError(format!(
                "admitted tensor constructor is missing its raw parser ports for {expression}"
            )));
        }
        let interface = self.inspect_interface(expression).ok_or_else(|| {
            StructureError::ParsingError(format!(
                "admitted tensor is missing its parser boundary for {expression}"
            ))
        })?;
        let slots = interface
            .logical_slots()
            .into_iter()
            .map(|slot| {
                let index = match slot.aind {
                    PartialIndex::Explicit(index) => convert_parser_index::<_, Aind>(index)?,
                    PartialIndex::Open(_) => next_index(),
                };
                Ok(slot.rep().slot(index))
            })
            .collect::<Result<Vec<_>, StructureError>>()?;
        Ok(OrderedStructure::new(slots))
    }

    /// A graph can use a caller's index grammar without admitting those
    /// identities into Idenso's independent PartialStructure representation.
    fn uses_custom_graph_indices<Aind: AbsInd + ParseableAind>(
        structure: &Canonicalized<OrderedStructure<LibraryRep, Aind>>,
    ) -> bool {
        std::any::TypeId::of::<Aind>() != std::any::TypeId::of::<AbstractIndex>()
            && structure
                .canonical()
                .external_structure_iter()
                .any(|slot| convert_parser_index::<_, AbstractIndex>(slot.aind()).is_err())
    }

    /// Raw graph parsing has its own tensor filter and logical slot layout.
    /// Once it has accepted a leaf, that layout establishes the raw incidences,
    /// including slots in an untagged constructor or projector.
    fn parsed_leaf<Aind: AbsInd + ParseableAind>(
        expression: AtomView<'_>,
        structure: &Canonicalized<OrderedStructure<LibraryRep, Aind>>,
        admission: Option<&Self>,
    ) -> Result<Self, TensorInferenceError> {
        if admission.is_some() {
            return Self::parser_leaf(expression, structure.canonical(), admission);
        }
        if Self::uses_custom_graph_indices(structure) {
            let mut proofs = Self {
                occurrences: OnceLock::from(Arc::default()),
                ..Self::default()
            };
            let raw = structure
                .canonical()
                .external_structure_iter()
                .map(|slot| slot.to_atom())
                .collect::<Vec<_>>();
            proofs
                .parser_leaf_interfaces()
                .lock()
                .expect("leaf interface cache is not poisoned")
                .insert(
                    expression.to_owned(),
                    (true, structure.layout().canonical_to_logical(&raw)),
                );
            proofs.publish_interfaces(expression, None, []);
            return Ok(proofs);
        }
        if !matches!(expression, AtomView::Fun(_)) {
            let raw = Self::raw_parser_regions(expression, StrictTensorFilter::ContainsReps)?;
            return Self::parser_leaf(expression, structure.canonical(), Some(&raw));
        }
        let supplied = Self::default();
        let ports = structure
            .canonical()
            .external_structure_iter()
            .map(|slot| slot.to_atom())
            .collect::<Vec<_>>();
        supplied.retain_leaf_ports(expression, structure.layout().canonical_to_logical(&ports))?;
        Self::parser_leaf(expression, structure.canonical(), Some(&supplied))
    }

    /// A raw parser may fold arithmetic into one opaque node. Admit its nested
    /// leaves under that parser's head filter, including legacy untagged slots.
    fn raw_parser_regions(
        expression: AtomView<'_>,
        filter: StrictTensorFilter,
    ) -> Result<Self, TensorInferenceError> {
        #[cfg(test)]
        network_metadata_tests::RAW_PARSER_ADMISSIONS.with(|count| count.set(count.get() + 1));
        let supplied = Self::default();
        let mut matcher = SlotMatcher::default();
        let mut pending = vec![expression];
        while let Some(node) = pending.pop() {
            if supplied.retained_interface(node).is_some() {
                continue;
            }
            let function = match node {
                AtomView::Fun(function) => function,
                AtomView::Add(sum) => {
                    pending.extend(sum.iter());
                    continue;
                }
                AtomView::Mul(product) => {
                    pending.extend(product.iter());
                    continue;
                }
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    pending.extend([base, exponent]);
                    continue;
                }
                _ => continue,
            };
            if !node.is_tensorial(filter) {
                supplied.retain_interface(node, PartialStructure::from_logical_slots([]));
                continue;
            }
            if composition::is_composite_head(function.get_symbol(), &SPENSO_TAG)
                || function.get_symbol() == SPENSO_TAG.dot
            {
                pending.extend(function.iter());
                continue;
            }
            let ports = if let Some(ports) =
                InterfaceInference::constructor_ports(function, &mut matcher, true)
            {
                ports
            } else {
                match OrderedStructure::<LibraryRep, AbstractIndex>::structure_from_atom(
                    node,
                    &mut matcher,
                ) {
                    Ok(structure) => structure.layout().canonical_to_logical(
                        &structure
                            .canonical()
                            .external_structure_iter()
                            .map(|slot| slot.to_atom())
                            .collect::<Vec<_>>(),
                    ),
                    Err(StructureError::EmptyStructure(_)) => Vec::new(),
                    Err(error) => return Err(TensorInferenceError::invalid(error.to_string())),
                }
            };
            supplied.retain_leaf_ports(node, ports)?;
        }
        Ok(supplied)
    }

    /// Attach the parser's established ports and reuse the admission owner's
    /// scoped facts. Only raw parser entry points admit unknown syntax here.
    fn parser_leaf<Aind: AbsInd + ParseableAind>(
        expression: AtomView<'_>,
        structure: &OrderedStructure<LibraryRep, Aind>,
        admission: Option<&Self>,
    ) -> Result<Self, TensorInferenceError> {
        let interface = Self::parser_boundary(structure)?;
        let Some(admission) = admission else {
            return SymbolicTensor::<PartialStructure>::from_observed_parts(
                expression.to_owned(),
                interface,
                LeafInference::ObserveEncoded,
            )
            .map(|tensor| tensor.proofs);
        };
        admission.admitted_parser_leaf(
            expression,
            structure,
            interface,
            admission.parser_history(expression),
        )
    }

    fn parser_history(&self, expression: AtomView<'_>) -> Option<&ExplicitIndexOccurrences> {
        self.occurrences
            .get()
            .filter(|_| self.owns_expression(expression))
            .map(Arc::as_ref)
    }

    fn admitted_parser_leaf<Aind: AbsInd + ParseableAind>(
        &self,
        expression: AtomView<'_>,
        structure: &OrderedStructure<LibraryRep, Aind>,
        interface: PartialStructure,
        inherited: Option<&ExplicitIndexOccurrences>,
    ) -> Result<Self, TensorInferenceError> {
        let (occurrences, local) =
            InterfaceInference::parser_occurrences(expression, self, Some(&interface))?;
        let mut proofs = local.unwrap_or_else(|| self.clone());
        if !self.owns_expression(expression) {
            // Product readiness belongs to the exact declared payload.
            proofs.normalized_product = OnceLock::new();
        }
        proofs.occurrences = OnceLock::from(Self::finish_parser_history(
            occurrences,
            structure,
            inherited,
        )?);
        proofs.publish_interfaces(expression, interface, []);
        Ok(proofs)
    }

    fn parser_boundary<Aind: AbsInd + ParseableAind>(
        structure: &OrderedStructure<LibraryRep, Aind>,
    ) -> Result<PartialStructure, TensorInferenceError> {
        let mut boundary = structure.clone();
        boundary.trace_out();
        Ok(PartialStructure::from_logical_slots(
            boundary
                .external_structure_iter()
                .map(|slot| {
                    convert_parser_index::<_, AbstractIndex>(slot.aind())
                        .map(|index| slot.rep().slot(PartialIndex::Explicit(index)))
                        .map_err(|error| TensorInferenceError::invalid(error.to_string()))
                })
                .collect::<Result<Vec<_>, _>>()?,
        ))
    }

    /// Count both written and stored ports, including exact-owner history that
    /// survived cancellation. All parser targets reject collisions eagerly.
    fn parser_occurrences<Aind: AbsInd + ParseableAind>(
        &self,
        expression: AtomView<'_>,
        structure: &OrderedStructure<LibraryRep, Aind>,
        inherited: Option<&ExplicitIndexOccurrences>,
    ) -> Result<Arc<ExplicitIndexOccurrences>, TensorInferenceError> {
        let (occurrences, _) = InterfaceInference::parser_occurrences(expression, self, None)?;
        Self::finish_parser_history(occurrences, structure, inherited)
    }

    fn finish_parser_history<Aind: AbsInd + ParseableAind>(
        mut occurrences: ExplicitIndexOccurrences,
        structure: &OrderedStructure<LibraryRep, Aind>,
        inherited: Option<&ExplicitIndexOccurrences>,
    ) -> Result<Arc<ExplicitIndexOccurrences>, TensorInferenceError> {
        if let Some(history) = inherited {
            occurrences.branch_max(history);
        }
        Self::with_parser_ports(Arc::new(occurrences), structure)
    }

    /// Merge the graph's storage axes with the literal's observed index uses.
    /// Reopening an occurrence reuses its checked history when all axes are covered.
    fn with_parser_ports<Aind: AbsInd + ParseableAind>(
        mut occurrences: Arc<ExplicitIndexOccurrences>,
        structure: &OrderedStructure<LibraryRep, Aind>,
    ) -> Result<Arc<ExplicitIndexOccurrences>, TensorInferenceError> {
        // The raw graph includes both the free boundary and any traced pairs.
        // Merge its counts once, preserving the literal's first written variance.
        let mut raw = ExplicitIndexOccurrences::default();
        for slot in structure.external_structure_iter() {
            let index = convert_parser_index::<_, AbstractIndex>(slot.aind())
                .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
            raw.add(slot.rep().slot(index), 1);
        }
        if !raw.is_bounded_by(&occurrences) {
            Arc::make_mut(&mut occurrences).branch_max(&raw);
        }
        occurrences.validate()?;
        Ok(occurrences)
    }
}

impl<S: PartialEq, E: PartialEq> PartialEq for SymbolicTensor<S, E> {
    fn eq(&self, other: &Self) -> bool {
        self.structure == other.structure
            && self.expression == other.expression
            && self.is_metric == other.is_metric
            && self.is_composite == other.is_composite
    }
}
impl<S: Eq, E: Eq> Eq for SymbolicTensor<S, E> {}

impl std::hash::Hash for SymbolicTensor<spenso::structure::partial::PartialStructure> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.hash_value(state);
    }
}

impl<S, E> SymbolicTensor<S, E> {
    pub fn expression(&self) -> &E {
        &self.expression
    }
    pub fn structure(&self) -> &S {
        &self.structure
    }
    pub fn is_metric(&self) -> bool {
        self.is_metric
    }
    pub fn is_composite(&self) -> bool {
        self.is_composite
    }

    /// Consume storage without exporting a reusable proof token.
    pub fn into_parts(self) -> (E, S) {
        (self.expression, self.structure)
    }
    pub fn into_expression(self) -> E {
        self.expression
    }

    /// Compare established reduction facts and retained cache identities.
    ///
    /// This is separate from tensor-value equality: an unchanged algebraic
    /// result can carry a new completion certificate. Callers reusing shared
    /// storage must check both. Cache identities are compared without walking
    /// their graphs or locking their contents.
    pub fn same_reduction_state(&self, other: &Self) -> bool {
        self.proofs.normalized_product.get() == other.proofs.normalized_product.get()
            && self.proofs.occurrences.get().map(Arc::as_ptr)
                == other.proofs.occurrences.get().map(Arc::as_ptr)
            && self.proofs.admitted_interfaces.get().map(Arc::as_ptr)
                == other.proofs.admitted_interfaces.get().map(Arc::as_ptr)
            && self
                .proofs
                .leaf_interfaces
                .same_source(&other.proofs.leaf_interfaces)
            && self.work.same_state(&other.work)
    }

    pub(crate) fn invalidate_proofs(&mut self) {
        self.replace_admission(TensorProofs::default());
    }

    /// A different admitted payload cannot reuse this tensor's reduction plan.
    pub(crate) fn replace_admission(&mut self, proofs: TensorProofs) {
        self.proofs = proofs;
        self.work = Default::default();
    }
}

impl<S, E: AtomCore + From<E::Output>> SymbolicTensor<S, E> {
    /// A sum can collapse to a product of indexed factors while its internal
    /// occurrences stay scoped to the original summands.
    pub(super) fn has_indexed_power_base(&self) -> Option<bool> {
        let history = self.proofs.occurrences.get()?;
        let open_factor = |factor: AtomView<'_>| {
            self.proofs
                .retained_interface(factor)
                .is_some_and(|interface| !interface.canonical().is_scalar())
        };
        Some(
            !history.is_empty()
                || match self.expression.as_atom_view() {
                    AtomView::Mul(product) => product.iter().any(open_factor),
                    AtomView::Pow(power) => open_factor(power.get_base_exp().0),
                    _ => false,
                },
        )
    }

    /// Protect an implicit scalar contraction when scalar arithmetic repeats it.
    /// Admitted callers provide their index metadata. Only legacy raw network
    /// storage without that metadata needs to discover its written scope.
    pub(super) fn scalar_power_expression(expression: E, indexed: Option<bool>) -> E {
        let indexed_product = matches!(
            expression.as_atom_view(),
            AtomView::Mul(_) | AtomView::Pow(_)
        ) && indexed.unwrap_or_else(|| {
            let mut indexed_product = false;
            expression.visitor(&mut |value| {
                if let AtomView::Fun(fun) = value {
                    let symbol = fun.get_symbol();
                    if symbol == SPENSO_TAG.bracket || symbol.is_scalar() {
                        return false;
                    }
                    indexed_product |=
                        symbol.has_attributes_of(SPENSO_TAG.rep_) && fun.get_nargs() == 2;
                }
                !indexed_product
            });
            indexed_product
        });
        // Keep implicit index sums atomic when scalar arithmetic takes powers.
        if indexed_product {
            expression
                .atom_to_output(bracket!(expression.as_atom_view()))
                .into()
        } else {
            expression
        }
    }
}

#[cfg(test)]
impl<Aind: AbsInd + ParseableAind> SymbolicTensor<OrderedStructure<LibraryRep, Aind>> {
    /// Observe a callback-sensitive result without changing its planned port order.
    pub(crate) fn validate_rewritten_interface(
        &self,
        expression: &Atom,
    ) -> Result<(), TensorInferenceError> {
        use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};

        if expression.as_view().is_zero() {
            return Ok(());
        }
        let slots = self
            .structure
            .external_structure()
            .into_iter()
            .map(|slot| {
                // Use the index's symbolic representation, not its Rust storage
                // type. The existing observer owns parsing and callback policy.
                let encoded = slot.aind().to_atom();
                let index = AbstractIndex::from_view(encoded.as_view()).map_err(|error| {
                    TensorInferenceError::Invalid(format!(
                        "cannot validate rewritten tensor index: {error}"
                    ))
                })?;
                if index.to_atom() != encoded {
                    return Err(TensorInferenceError::Invalid(
                        "rewritten tensor index does not have a lossless shared representation"
                            .into(),
                    ));
                }
                Ok(slot.rep().slot(PartialIndex::Explicit(index)))
            })
            .collect::<Result<Vec<_>, TensorInferenceError>>()?;
        SymbolicTensor::<PartialStructure>::validate_encoded_interface(
            expression,
            &PartialStructure::from_logical_slots(slots),
        )
    }
}

impl<S, E> Ref for SymbolicTensor<S, E> {
    type Ref<'a>
        = &'a SymbolicTensor<S, E>
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        self
    }
}

impl<Aind: AbsInd> spenso::network::FastTensorSum
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
}
impl<Aind: AbsInd> spenso::network::FastTensorSumContractible<Atom>
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
}
impl<Aind: AbsInd> spenso::network::TensorCommonFactor<Atom>
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
}

impl<Aind: AbsInd> ScalarTensor for SymbolicTensor<OrderedStructure<LibraryRep, Aind>> {
    fn scalar_power(self, exponent: i8, power: impl FnOnce(Atom) -> Atom) -> Self {
        assert!(
            self.structure.is_scalar(),
            "scalar power requires a scalar tensor"
        );
        let occurrences = self
            .proofs
            .occurrences
            .get()
            .expect("published scalar tensor is missing construction history");
        let expression = if exponent == 1 {
            self.expression.clone()
        } else {
            Self::scalar_power_expression(self.expression.clone(), self.has_indexed_power_base())
        };
        let expression = power(expression);
        let mut proofs = self.proofs.clone();
        if exponent != 1 {
            // Every power of a scalar binds its internal contractions. Reusing
            // the carrier must not expose those local identities to a product.
            proofs.occurrences = OnceLock::from(Arc::new(
                occurrences.scoped(&PartialStructure::from_logical_slots([])),
            ));
        }
        proofs.publish_interfaces(
            expression.as_view(),
            PartialStructure::from_logical_slots([]),
            [],
        );
        Self {
            expression,
            work: Default::default(),
            proofs,
            is_composite: true,
            is_metric: false,
            structure: self.structure,
        }
    }

    /// # Panics
    /// Panics if the supplied scalar has invalid explicit-index multiplicities.
    fn new_scalar(scalar: Self::Scalar) -> Self {
        if InterfaceInference::has_structured_syntax(scalar.as_view()) {
            let raw =
                TensorProofs::raw_parser_regions(scalar.as_view(), StrictTensorFilter::Tagged)
                    .expect("scalar tensor intake requires valid parser ports");
            let proofs = TensorProofs::parser_leaf(
                scalar.as_view(),
                &OrderedStructure::<LibraryRep, AbstractIndex>::scalar_structure(),
                Some(&raw),
            )
            .expect("scalar tensor intake requires valid scalar index scopes");
            return Self {
                expression: scalar,
                structure: OrderedStructure::scalar_structure(),
                work: Default::default(),
                proofs,
                is_metric: false,
                is_composite: true,
            };
        }
        // This is fresh scalar intake, not a transfer of an existing carrier.
        // Structured scalar payloads can retain internal pairs even when their
        // free interface is empty; do not replace that history by an empty mask.
        let occurrences =
            ExplicitIndexOccurrences::from_atom(scalar.as_view(), &mut SlotMatcher::default());
        occurrences
            .validate()
            .expect("scalar tensor producer requires valid index multiplicities");
        let mut proofs = TensorProofs {
            occurrences: OnceLock::from(Arc::new(occurrences)),
            ..TensorProofs::default()
        };
        proofs.publish_interfaces(
            scalar.as_view(),
            PartialStructure::from_logical_slots([]),
            [],
        );
        SymbolicTensor {
            work: Default::default(),
            proofs,
            structure: OrderedStructure::scalar_structure(),
            is_metric: false,
            is_composite: false,
            expression: scalar,
        }
    }
}

impl<Aind: AbsInd> ScalarStructure for SymbolicTensor<OrderedStructure<LibraryRep, Aind>> {
    fn scalar_structure() -> Self {
        Self::new_scalar(Atom::num(1))
    }
}

impl<Aind: AbsInd + DummyAind + ParseableAind> TensorIdentity
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
    type Id = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>;
    type IdSlot = LibrarySlot<Aind>;

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        Self::from_named(&NamedStructure::<Symbol, (), LibraryRep, Aind>::id(i, j)).unwrap()
    }
}

impl<Aind: AbsInd + DummyAind + ParseableAind> ApplyPendingIndexPermutation
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
    type Output = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>;

    fn apply_pending_index_permutation(self, pending: &PendingIndexPermutation) -> Self::Output {
        if pending.is_identity() {
            return self;
        }

        let indices = pending
            .apply_slice_inverse(&self.structure.external_structure())
            .into_iter()
            .map(|slot| slot.aind())
            .collect::<Vec<_>>();
        // Symbolic coordinate order lives in the written tensor arguments.
        // Relabel those occurrences together with their retained contracts.
        self.reindex_storage(&indices)
            .expect("checked index permutation must preserve tensor ports")
            .apply()
    }
}

impl<Aind: AbsInd + ParseableAind, E: AtomCore + From<E::Output>>
    SymbolicTensor<OrderedStructure<LibraryRep, Aind>, E>
{
    /// Expose storage-axis occurrences to the existing port transformer. A
    /// self-loop has two storage axes even though its public boundary is closed.
    fn storage_port_tensor(
        &self,
    ) -> Result<SymbolicTensor<PartialStructure>, TensorInferenceError> {
        let slots = self
            .structure
            .external_structure_iter()
            .map(|slot| {
                convert_parser_index::<_, AbstractIndex>(slot.aind())
                    .map(|index| slot.rep().slot(PartialIndex::Explicit(index)))
                    .map_err(|error| TensorInferenceError::invalid(error.to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let structure = PartialStructure::from_logical_slots(slots);
        let expression = self.expression.as_atom_view().to_owned();
        let mut proofs = self.proofs.clone();
        proofs.publish_interfaces(expression.as_view(), structure.clone(), []);
        let source = SymbolicTensor {
            structure,
            expression,
            work: Default::default(),
            proofs,
            is_metric: self.is_metric,
            is_composite: self.is_composite,
        };
        source.explicit_index_occurrences()?;
        Ok(source)
    }
}

impl<Aind: AbsInd + ParseableAind, E: AtomCore + From<E::Output>> TensorStructure
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>, E>
{
    // type R = <T::Structure as TensorStructure>::R;
    type Indexed = SymbolicTensor<OrderedStructure<LibraryRep, Aind>, E>;
    type Slot = LibrarySlot<Aind>;

    fn reindex_storage(self, indices: &[Aind]) -> Result<Reindexed<Self::Indexed>, StructureError> {
        let target = self.structure.clone().reindex_storage(indices)?.apply();
        let source = self
            .storage_port_tensor()
            .map_err(|error| StructureError::ParsingError(error.to_string()))?;
        let replacements = indices
            .iter()
            .enumerate()
            .map(|(position, index)| {
                convert_parser_index::<_, AbstractIndex>(*index).map(|index| (position, index))
            })
            .collect::<Result<HashMap<_, _>, _>>()
            .map_err(|error| StructureError::ParsingError(error.to_string()))?;
        let mapped = source
            .reindex_interface_ports(&replacements)
            .map_err(|error| StructureError::ParsingError(error.to_string()))?;
        let interface = InterfaceInference::merge_explicit_interface_sequence(&[mapped.structure])
            .map_err(|error| StructureError::ParsingError(error.to_string()))?;
        let mut proofs = mapped.proofs;
        proofs.publish_interfaces(mapped.expression.as_view(), interface, []);
        let expression = self.expression.atom_to_output(mapped.expression).into();
        let target_indices = target.external_indices_iter().collect::<Vec<_>>();
        // The expression already carries the coordinate permutation. Return
        // an identity transition so applying Reindexed does not repeat it.
        Ok(target
            .reindex_storage(&target_indices)?
            .map_target(|structure| Self {
                structure,
                expression,
                work: Default::default(),
                proofs,
                is_metric: self.is_metric,
                is_composite: self.is_composite,
            }))
    }

    fn dual(self) -> Self {
        let source = self
            .storage_port_tensor()
            .expect("dual requires an admitted symbolic tensor");
        let ports = source
            .structure
            .logical_slots()
            .into_iter()
            .collect::<HashSet<_>>();
        let map_interface = |interface: &PartialStructure| {
            Ok(PartialStructure::from_logical_slots(
                interface.logical_slots().into_iter().map(|slot| {
                    if ports.contains(&slot) {
                        slot.rep().dual().slot(slot.aind)
                    } else {
                        slot
                    }
                }),
            ))
        };
        let mut matcher = SlotMatcher::default();
        let (expression, mut proofs) = source
            .map_admitted_indices(
                |value| {
                    Ok(matcher
                        .parse::<LibraryRep, AbstractIndex>(value)
                        .ok()
                        .and_then(|slot| {
                            ports
                                .contains(&slot.rep().slot(PartialIndex::Explicit(slot.aind())))
                                .then(|| {
                                    slot.rep()
                                        .dual()
                                        .slot::<AbstractIndex, _>(slot.aind())
                                        .to_atom()
                                })
                        }))
                },
                map_interface,
            )
            .expect("dual must preserve admitted tensor syntax");
        let occurrences = source
            .explicit_index_occurrences()
            .unwrap()
            .mapped(|slot| {
                Ok::<_, TensorInferenceError>(
                    if ports.contains(&slot.rep().slot(PartialIndex::Explicit(slot.aind()))) {
                        slot.rep().dual().slot(slot.aind())
                    } else {
                        slot
                    },
                )
            })
            .expect("dual preserves index identities");
        occurrences
            .validate()
            .expect("dual preserves valid index multiplicities");
        proofs.occurrences = OnceLock::from(Arc::new(occurrences));

        let interface = InterfaceInference::merge_explicit_interface_sequence(&[map_interface(
            &source.structure,
        )
        .unwrap()])
        .expect("dual preserves compatible port pairs");
        TensorProofs::assert_normalization_preserves_interface(
            "dual",
            &[&source.expression],
            &expression,
            Ok(&interface),
            None,
        );
        proofs.publish_interfaces(expression.as_view(), interface, []);
        Self {
            expression: self.expression.atom_to_output(expression).into(),
            structure: self.structure.dual(),
            work: Default::default(),
            proofs,
            is_metric: self.is_metric,
            is_composite: self.is_composite,
        }
    }

    delegate! {
        to self.structure() {
            fn is_fully_self_dual(&self)-> bool;
            fn external_reps_iter(&self)-> impl Iterator<Item = Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn external_indices_iter(&self)-> impl Iterator<Item = Aind>;
            fn external_dims_iter(&self)-> impl Iterator<Item = Dimension>;
            fn external_structure_iter(&self)-> impl Iterator<Item = Self::Slot>;
            fn get_slot(&self,i:impl Into<SlotIndex>)-> Option<Self::Slot>;
            fn get_rep(&self,i:impl Into<SlotIndex>)-> Option<Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn get_dim(&self,i:impl Into<SlotIndex>)-> Option<Dimension>;
            fn get_aind(&self,i:impl Into<SlotIndex>)-> Option<Aind>;
            fn order(&self)-> usize;
        }
    }
}

impl<S: TensorStructure, E: AtomCore + From<E::Output>> HasStructure for SymbolicTensor<S, E> {
    const PRESERVE_SCALAR_METADATA: bool = true;

    type Structure = S;
    type Scalar = E;
    type ScalarRef<'a>
        = &'a E
    where
        Self: 'a;
    type Store<O>
        = SymbolicTensor<O, E>
    where
        O: TensorStructure;

    fn map_structure<O: TensorStructure>(
        self,
        f: impl FnOnce(Self::Structure) -> O,
    ) -> Self::Store<O> {
        SymbolicTensor {
            work: Default::default(),
            proofs: Default::default(),
            structure: f(self.structure),
            expression: self.expression,
            is_metric: self.is_metric,
            is_composite: self.is_composite,
        }
    }

    fn map_structure_result<O: TensorStructure, Er>(
        self,
        f: impl FnOnce(Self::Structure) -> eyre::Result<O, Er>,
    ) -> std::result::Result<Self::Store<O>, Er> {
        Ok(SymbolicTensor {
            work: Default::default(),
            proofs: Default::default(),
            structure: f(self.structure)?,
            expression: self.expression,
            is_metric: self.is_metric,
            is_composite: self.is_composite,
        })
    }

    fn structure(&self) -> &Self::Structure {
        &self.structure
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        self.invalidate_proofs();
        &mut self.structure
    }

    fn map_same_structure(self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        SymbolicTensor {
            work: Default::default(),
            proofs: Default::default(),
            structure: f(self.structure),
            is_metric: self.is_metric,
            is_composite: self.is_composite,
            expression: self.expression,
        }
    }

    fn scalar(self) -> Option<Self::Scalar> {
        if self.structure.is_scalar() {
            Some(self.expression)
        } else {
            None
        }
    }

    fn scalar_power_base(self) -> Option<Self::Scalar> {
        let indexed = self.has_indexed_power_base();
        self.scalar()
            .map(|expression| Self::scalar_power_expression(expression, indexed))
    }

    fn scalar_ref(&self) -> Option<&Self::Scalar> {
        if self.structure.is_scalar() {
            Some(&self.expression)
        } else {
            None
        }
    }
}

impl<Aind: AbsInd + ParseableAind> IteratableTensor
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
    type Data<'a>
        = AtomView<'a>
    where
        Self: 'a;

    fn iter_flat(&self) -> impl Iterator<Item = (FlatIndex, Self::Data<'_>)> {
        std::iter::once((FlatIndex::from(0usize), self.expression.as_view()))
    }

    fn iter_expanded(&self) -> impl Iterator<Item = (ExpandedIndex, Self::Data<'_>)> {
        std::iter::once((
            ExpandedIndex {
                indices: Vec::new(),
            },
            self.expression.as_view(),
        ))
    }
}

impl<Aind: AbsInd> StructureContract for SymbolicTensor<OrderedStructure<LibraryRep, Aind>> {
    fn merge(
        &self,
        other: &Self,
    ) -> Result<(Self, SubSet<SlotIndex>, SubSet<SlotIndex>, MergeInfo), StructureError> {
        let expression = &other.expression * &self.expression;
        let (structure, pos_self, pos_other, mergeinfo) = self.structure.merge(&other.structure)?;
        let left = self.proofs.retained_interface(self.expression.as_view());
        let right = other.proofs.retained_interface(other.expression.as_view());
        let left_history = self
            .proofs
            .occurrences
            .get()
            .expect("published tensor is missing construction history");
        let right_history = other
            .proofs
            .occurrences
            .get()
            .expect("published tensor is missing construction history");
        // Graph edges establish the contracted pairings. An index bound by
        // an earlier contraction must not attach to an independent later edge
        // that uses the same label, including powers written as repeated leaves.
        let (interface, mut occurrences) = if let (Some(left), Some(right)) = (left, right) {
            let mut active = left_history.scoped(&left);
            active.append(&right_history.scoped(&right));
            let interface = InterfaceInference::merge_explicit_interface_sequence(&[left, right])
                .map_err(|error| StructureError::ParsingError(error.to_string()))?;
            (Some(interface), active)
        } else {
            (None, ExplicitIndexOccurrences::default())
        };
        occurrences.branch_max(left_history);
        occurrences.branch_max(right_history);
        occurrences
            .validate()
            .map_err(|error| StructureError::ParsingError(error.to_string()))?;
        let mut proofs = TensorProofs {
            occurrences: OnceLock::from(Arc::new(occurrences)),
            ..TensorProofs::default()
        };
        proofs.publish_interfaces(
            expression.as_view(),
            interface,
            [&self.proofs, &other.proofs],
        );

        Ok((
            Self {
                work: Default::default(),
                proofs,
                structure,
                is_composite: true,
                is_metric: false,
                expression,
            },
            pos_self,
            pos_other,
            mergeinfo,
        ))
    }

    fn trace_out(&mut self) {
        self.work = Default::default();
        self.structure.trace_out();
    }

    fn trace(&mut self, i: usize, j: usize) {
        self.work = Default::default();
        self.structure.trace(i, j);
    }
}

impl<Aind: AbsInd> Trace for SymbolicTensor<OrderedStructure<LibraryRep, Aind>> {
    fn internal_contract(&self) -> Self {
        let mut traced = self.clone();
        traced.trace_out();
        traced.is_composite = true;
        traced.is_metric = false;
        traced
    }
}

// impl<Const> Shadowable<Const> for SymbolicTensor {}

#[allow(dead_code)]
impl<Aind: AbsInd> SymbolicTensor<OrderedStructure<LibraryRep, Aind>> {
    pub fn from_named<N>(structure: &N) -> Option<Self>
    where
        Aind: ParseableAind,
        N: ToSymbolic + HasName + TensorStructure<Slot = LibrarySlot<Aind>>,
        N::Name: IntoSymbol + Clone,
        N::Args: IntoArgs,
    {
        let canonicalized = Canonicalized::from(structure.external_structure());
        let name = structure.name()?.ref_into_symbol();
        let expression = structure.to_symbolic(None)?;
        let proofs = Self::constructor_proofs(
            &expression,
            canonicalized.canonical(),
            name,
            structure
                .external_structure_iter()
                .map(|slot| slot.to_atom())
                .collect(),
        )
        .ok()?;
        Some(SymbolicTensor {
            work: Default::default(),
            proofs,
            expression,
            is_metric: name == ETS.metric,
            is_composite: false,
            structure: canonicalized.into_canonical(),
        })
    }

    pub fn from_canonicalized<N>(structure: &Canonicalized<N>) -> Option<Self>
    where
        Aind: ParseableAind,
        N: ToSymbolic + HasName + TensorStructure<Slot = LibrarySlot<Aind>>,
        N::Name: IntoSymbol + Clone,
        N::Args: IntoArgs,
    {
        let canonical = structure.canonical();
        let canonicalized = Canonicalized::from(canonical.external_structure());
        let name = canonical.name()?.ref_into_symbol();
        let expression = structure.to_symbolic(None)?;
        let logical_ports = structure.layout().canonical_to_logical(
            &canonical
                .external_structure_iter()
                .map(|slot| slot.to_atom())
                .collect::<Vec<_>>(),
        );
        let proofs =
            Self::constructor_proofs(&expression, canonicalized.canonical(), name, logical_ports)
                .ok()?;
        Some(SymbolicTensor {
            work: Default::default(),
            proofs,
            expression,
            is_composite: false,
            is_metric: name == ETS.metric,
            structure: canonicalized.into_canonical(),
        })
    }

    fn constructor_proofs(
        expression: &Atom,
        structure: &OrderedStructure<LibraryRep, Aind>,
        head: Symbol,
        mut ports: Vec<Atom>,
    ) -> Result<TensorProofs, TensorInferenceError>
    where
        Aind: ParseableAind,
    {
        if !InterfaceInference::intrinsic_normalization_head(head) {
            // A user normalizer may rewrite the leaf's internal incidences.
            // Its free-interface contract alone does not certify raw graph ports.
            return TensorProofs::parser_leaf(expression.as_view(), structure, None);
        }
        let proofs = TensorProofs::parser_leaf(
            expression.as_view(),
            structure,
            Some(&TensorProofs::default()),
        )?;
        if let AtomView::Fun(function) = expression.as_view() {
            let mut written = Vec::with_capacity(ports.len());
            if let Some(first_port) = function.get_nargs().checked_sub(ports.len()) {
                for argument in function.iter().skip(first_port) {
                    let Some(position) = ports.iter().position(|port| port.as_view() == argument)
                    else {
                        break;
                    };
                    written.push(ports.remove(position));
                }
            }
            if function.get_symbol() == head && ports.is_empty() {
                // Intrinsic symmetries can reorder their constructor arguments;
                // match known typed ports without decoding the emitted syntax.
                proofs.retain_leaf_ports(expression.as_view(), written)?;
            } else if let Some(leaves) = proofs.leaf_interfaces.get() {
                leaves
                    .lock()
                    .expect("leaf interface cache is not poisoned")
                    .remove(expression.as_view().get_data());
            }
        }
        Ok(proofs)
    }

    pub fn to_named(&self) -> NamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind> {
        NamedStructure {
            structure: self.structure.clone(),
            global_name: self.name(),
            additional_args: self.args(),
        }
    }

    pub fn empty(expression: Atom) -> Self {
        SymbolicTensor {
            work: Default::default(),
            proofs: Default::default(),
            structure: OrderedStructure::empty(),
            is_composite: false,
            is_metric: false,
            expression,
        }
    }

    // #[must_use]
    pub fn get_atom(&self) -> &Atom {
        &self.expression
    }

    // pub fn to_mixed(self) -> MixedTensor {
    //     self.to_named().to_shell().to_explicit().unwrap()
    // }
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn to_network(
        &self,
        library: &TensorLibrary<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>,
    ) -> Result<
        Network<
            NetworkStore<MixedTensor<f64, ShadowedStructure<AbstractIndex>>, Atom>,
            ExplicitKey<AbstractIndex>,
            Symbol,
            AbstractIndex,
        >,
        TensorNetworkError<ExplicitKey<AbstractIndex>, Symbol>,
    > {
        let settings =
            ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
        type ComponentNet = Network<
            NetworkStore<MixedTensor<f64, ShadowedStructure<AbstractIndex>>, Atom>,
            ExplicitKey<AbstractIndex>,
            Symbol,
            AbstractIndex,
        >;
        if self.proofs.is_admitted() {
            ComponentNet::try_from_admitted_view_with_context(
                self.expression.as_view(),
                library,
                &settings,
                self.proofs.leaf_interfaces.read_cache().unwrap_or_default(),
                None,
            )
        } else {
            ComponentNet::try_from_view(self.expression.as_view(), library, &settings)
        }
    }
}

// impl TryFrom<Atom> for SymbolicTensor {
//     type Error = String;
//     fn try_from(value: Atom) -> Result<Self, Self::Error> {
//         let structure = value
//             .as_view()
//             .try_into()
//             .unwrap_or(OrderedStructure::empty());

//         Ok(SymbolicTensor {
//             structure,
//             expression: value,
//         })
//     }
// }

impl<S> HasName for SymbolicTensor<S> {
    type Name = Symbol;
    type Args = Vec<Atom>;

    fn name(&self) -> Option<Self::Name> {
        if let AtomView::Fun(f) = self.expression.as_view() {
            Some(f.get_symbol())
        } else {
            None
        }
    }

    fn set_name(&mut self, _name: Self::Name) {
        unimplemented!("Cannot set name of a symbolic tensor")
    }

    fn args(&self) -> Option<Self::Args> {
        let mut args = vec![];
        match self.expression.as_view() {
            AtomView::Fun(f) => {
                for arg in f.iter() {
                    if let AtomView::Fun(f) = arg {
                        if f.get_symbol() != AIND_SYMBOLS.aind {
                            args.push(arg.to_owned());
                        }
                    } else {
                        args.push(arg.to_owned());
                    }
                }
                Some(args)
            }
            _ => None,
        }
    }
}

/// Symbolic contraction of two symbolic tensors is just a multiplication of the atoms.
///
impl<Aind: AbsInd> Contract<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
    type LCM = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>;
    fn contract(
        &self,
        other: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    ) -> Result<Self::LCM, ContractionError> {
        let (out, _, _, _) = self.merge(other)?;
        Ok(out)
    }
}

pub type SymbolicNet<Aind, E = Atom> = Network<
    NetworkStore<SymbolicTensor<OrderedStructure<LibraryRep, Aind>, E>, E>,
    DummyKey,
    Symbol,
    Aind,
>;

// pub type HepTensor<Aind> = MixedTensor<f64, ShadowedStructure<Aind>>;

// pub type ParamNet<Aind> =
//     Network<NetworkStore<ParamTensor<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom>, DummyKey, Symbol, Aind>;
pub trait SymbolicNetExt<Aind: AbsInd + DummyAind + ParseableAind + 'static> {
    fn snapshot_dot(&self) -> String;
    /// Return whether signed-zero pruning changed the existing graph.
    fn remove_antisymmetric_zero_terms(&mut self) -> bool;
    fn simple_execute<CStrat>(self) -> Result<Atom, NetworkToolingError>
    where
        SymbolicTensor<OrderedStructure<LibraryRep, Aind>>: Contract<
                SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
                CStrat,
                LCM = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
            >;
}

impl<Aind: AbsInd + DummyAind + ParseableAind + 'static> SymbolicNetExt<Aind>
    for SymbolicNet<Aind>
{
    fn remove_antisymmetric_zero_terms(&mut self) -> bool {
        canonicalize::remove_antisymmetric_zero_terms(self)
    }

    fn snapshot_dot(&self) -> String {
        self.dot_display_impl(
            |a| a.to_bare_ordered_string(),
            |_| None,
            |a| a.expression.to_bare_ordered_string(),
            |a| a.to_string(),
        )
    }

    fn simple_execute<CStrat>(mut self) -> Result<Atom, NetworkToolingError>
    where
        SymbolicTensor<OrderedStructure<LibraryRep, Aind>>: Contract<
                SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
                CStrat,
                LCM = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
            >,
    {
        let lib = DummyLibrary::<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>::new();

        self.execute::<Sequential, SmallestDegree<CStrat>, _, _, _>(&lib, &Wrap {})
            .map_err(|error| NetworkToolingError::Execute {
                reason: error.to_string(),
            })?;

        Ok(
            match self
                .result_tensor(&lib)
                .map_err(|error| NetworkToolingError::Result {
                    reason: error.to_string(),
                })? {
                ExecutionResult::One => Atom::num(1),
                ExecutionResult::Zero => Atom::Zero,
                ExecutionResult::Val(tensor) => tensor.expression.clone(),
            },
        )
    }
}

pub trait SymbolicNetParse {
    #[allow(clippy::result_large_err)]
    fn parse_to_symbolic_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>>;
}

impl SymbolicNetParse for Atom {
    fn parse_to_symbolic_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        self.as_view().parse_to_symbolic_net::<Aind>(settings)
    }
}

impl SymbolicNetParse for AtomView<'_> {
    fn parse_to_symbolic_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        let lib = DummyLibrary::<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>::new();
        let settings = settings
            .clone()
            .with_strict_tensor_filter(StrictTensorFilter::ContainsReps);

        SymbolicNet::<Aind>::try_from_view::<OrderedStructure<LibraryRep, Aind>, _>(
            *self, &lib, &settings,
        )
    }
}

pub mod contract;
pub use contract::ContractSettings;
pub use simplification::{AlgebraContraction, AlgebraSettings, ReductionStatus};

impl<Aind: AbsInd + ParseableAind> Concretize<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
    for ShadowedStructure<Aind>
{
    fn concretize(self) -> eyre::Result<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>> {
        let name = self
            .name()
            .ok_or_else(|| eyre!("symbolic tensor requires a name"))?;
        let expression = self.to_symbolic(None).expect("the constructor has a name");
        let proofs = SymbolicTensor::constructor_proofs(
            &expression,
            &self.structure,
            name,
            self.external_structure_iter()
                .map(|slot| slot.to_atom())
                .collect(),
        )?;
        Ok(SymbolicTensor {
            work: Default::default(),
            proofs,
            expression,
            is_composite: false,
            is_metric: name == ETS.metric,
            structure: self.structure,
        })
    }

    fn concretize_logical(
        self,
        layout: &CanonicalLayout,
    ) -> eyre::Result<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>> {
        let name = self
            .name()
            .ok_or_else(|| eyre!("symbolic tensor requires a name"))?;
        let logical_ports = layout
            .canonical_to_logical(&self.external_structure())
            .into_iter()
            .map(|slot| slot.to_atom())
            .collect::<Vec<_>>();
        let expression = FunctionBuilder::new(name)
            .add_args(self.args().unwrap_or_default())
            .add_args(&logical_ports)
            .finish();
        let proofs =
            SymbolicTensor::constructor_proofs(&expression, &self.structure, name, logical_ports)?;
        Ok(SymbolicTensor {
            work: Default::default(),
            proofs,
            expression,
            is_composite: false,
            is_metric: name == ETS.metric,
            structure: self.structure,
        })
    }
}

impl<Aind: AbsInd + ParseableAind> Concretize<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
    for TensorShell<ShadowedStructure<Aind>>
{
    fn concretize(self) -> eyre::Result<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>> {
        self.structure.concretize()
    }

    fn concretize_logical(
        self,
        layout: &CanonicalLayout,
    ) -> eyre::Result<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>> {
        self.structure.concretize_logical(layout)
    }
}

impl<Aind: AbsInd + ParseableAind> Concretize<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
    for TensorShell<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
{
    fn concretize(self) -> eyre::Result<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>> {
        Ok(self.structure)
    }

    fn concretize_logical(
        self,
        layout: &CanonicalLayout,
    ) -> eyre::Result<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>> {
        let positions = (0..self.structure.order()).collect::<Vec<_>>();
        if layout.logical_to_canonical(&positions) != positions {
            eyre::bail!("cannot apply a logical layout to an already concrete symbolic tensor");
        }
        Ok(self.structure)
    }
}

impl<S: TensorInterface> std::ops::Neg for SymbolicTensor<S> {
    type Output = SymbolicTensor<S>;

    fn neg(self) -> Self::Output {
        let interface = self.proofs.retained_interface(self.expression.as_view());
        let expression = -self.expression;
        let mut proofs = self.proofs.clone();
        if let Some(interface) = interface {
            proofs.publish_interfaces(expression.as_view(), interface, []);
        }
        Self {
            work: Default::default(),
            proofs,
            expression,
            is_composite: true,
            is_metric: self.is_metric,
            structure: self.structure.into_algebra_result(),
        }
    }
}

impl<Aind: AbsInd> AddAssign<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
    fn add_assign(&mut self, rhs: SymbolicTensor<OrderedStructure<LibraryRep, Aind>>) {
        *self += &rhs;
    }
}

impl<Aind: AbsInd> AddAssign<&SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
    fn add_assign(&mut self, rhs: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>) {
        debug_assert_eq!(self.structure, rhs.structure);
        let interface = self.proofs.retained_interface(self.expression.as_view());
        let rhs_interface = rhs.proofs.retained_interface(rhs.expression.as_view());
        let mut occurrences = ExplicitIndexOccurrences::default();
        if let (Some(left), Some(right)) = (&interface, &rhs_interface) {
            occurrences = self
                .proofs
                .occurrences
                .get()
                .expect("published tensor is missing construction history")
                .scoped(left);
            occurrences.branch_max(
                &rhs.proofs
                    .occurrences
                    .get()
                    .expect("published tensor is missing construction history")
                    .scoped(right),
            );
        }
        let mut proofs = TensorProofs {
            occurrences: OnceLock::from(Arc::new(occurrences)),
            ..TensorProofs::default()
        };
        self.expression += &rhs.expression;
        proofs.publish_interfaces(
            self.expression.as_view(),
            interface,
            [&self.proofs, &rhs.proofs],
        );
        self.replace_admission(proofs);
        self.is_composite = true;
        self.is_metric = false;
    }
}

impl<Aind: AbsInd> ScalarMul<Atom> for SymbolicTensor<OrderedStructure<LibraryRep, Aind>> {
    type Output = Self;

    fn scalar_mul(&self, rhs: &Atom) -> Option<Self::Output> {
        self.merge(&Self::new_scalar(rhs.clone()))
            .ok()
            .map(|(tensor, ..)| tensor)
    }
}

impl ScalarMul<Atom> for SymbolicTensor<PartialStructure> {
    type Output = Self;

    /// # Panics
    ///
    /// Panics on invalid inputs; use [`Self::try_scalar_mul`] for typed errors.
    fn scalar_mul(&self, rhs: &Atom) -> Option<Self> {
        Some(self.try_scalar_mul(rhs).unwrap_or_else(|error| {
            panic!("symbolic tensor scalar multiplication failed: {error}")
        }))
    }
}

impl ScalarMul<Atom> for SymbolicTensor<IndexedStructure> {
    type Output = Self;

    /// # Panics
    ///
    /// Panics on invalid inputs; use [`Self::try_scalar_mul`] for typed errors.
    fn scalar_mul(&self, rhs: &Atom) -> Option<Self> {
        Some(self.try_scalar_mul(rhs).unwrap_or_else(|error| {
            panic!("indexed symbolic tensor scalar multiplication failed: {error}")
        }))
    }
}

impl<S: std::fmt::Display> std::fmt::Display for SymbolicTensor<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{},{}",
            self.structure,
            self.expression
                .printer(SpensoPrintSettings::compact().nice_symbolica())
        )
    }
}

#[cfg(test)]
#[test]
fn structure_mapping_keeps_expression_and_typed_zero() {
    use spenso::structure::representation::{Minkowski, RepName};
    let structure = OrderedStructure::<LibraryRep, AbstractIndex>::new(vec![
        Minkowski::default()
            .new_rep(Dimension::Concrete(4))
            .to_lib()
            .slot(AbstractIndex::Normal(431)),
    ])
    .into_canonical();
    let value = SymbolicTensor {
        work: Default::default(),
        proofs: Default::default(),
        expression: Atom::Zero,
        structure,
        is_metric: false,
        is_composite: true,
    };
    let name = symbolica::symbol!("mapped_symbolic_tensor");
    let mapped = value.map_structure(|structure| NamedStructure {
        structure,
        global_name: Some(name),
        additional_args: Some(Vec::<Atom>::new()),
    });
    assert!(mapped.expression.as_view().is_zero());
    assert_eq!(mapped.structure.order(), 1);
    assert_eq!(mapped.structure.global_name, Some(name));
    assert!(mapped.scalar_ref().is_none());
    let restored = mapped
        .map_structure_result(|named| Ok::<_, ()>(named.structure))
        .unwrap();
    assert!(restored.expression.as_view().is_zero());
    assert!(restored.is_composite);
    assert!(!restored.is_metric);
}
