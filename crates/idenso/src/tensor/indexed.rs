//! Fully indexed public boundaries backed by the existing tensor algebra.
//!
//! Owned engine projection moves the payload, layout and proof state. Borrowed
//! checked-result publication and dummy reservation use the original expression
//! and proof cells. Reduction owns one source copy and moves it through the shared
//! planner. Remaining borrowed adapters still clone the full owned Atom before
//! invoking their existing engine operations; index metadata is also O(rank).
//! No conversion re-infers an admitted payload or discards completion facts.
use super::composition::TensorCompositionError;
use super::inference::{BorrowedTensorParts, TensorInferenceError};
use super::simplification::{ReductionRequest, ReductionStatus};
use super::{AlgebraSettings, ContractSettings, SymbolicTensor};
use crate::{AdjointError, CanonicalizationError};
#[cfg(feature = "bincode")]
use bincode::de::Decoder;
#[cfg(feature = "bincode")]
use bincode::enc::Encoder;
#[cfg(feature = "bincode")]
use bincode::error::{DecodeError, EncodeError};
#[cfg(feature = "bincode")]
use bincode::{Decode, Encode};
use spenso::network::library::FunctionLibrary;
use spenso::network::parsing::{
    ParseSettings, ParseState, StructureFromAtom, TensorFromExpression, TensorLibraryFor,
};
use spenso::network::store::NetworkStore;
use spenso::network::{Network, TensorNetworkError};
use spenso::shadowing::TensorCollectFilter;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::dimension::Dimension;
use spenso::structure::partial::PartialStructure;
use spenso::structure::representation::{LibraryRep, LibrarySlot};
use spenso::structure::slot::{AbsInd, DummyAind, ParseableAind};
use spenso::structure::{
    Canonicalized, HasStructure, IndexedStructure, OrderedStructure, ScalarStructure,
    TensorStructure,
};
use spenso::tensors::data::SparseTensor;
use std::collections::HashMap;
use std::hash::Hasher;
#[cfg(feature = "bincode")]
use std::result::Result as StdResult;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::id::ConditionResult;
#[cfg(feature = "bincode")]
use symbolica::state::HasStateMap;

type IndexedTensor = SymbolicTensor<IndexedStructure>;
type PartialTensor = SymbolicTensor<PartialStructure>;
type Result<T> = std::result::Result<T, TensorInferenceError>;

impl From<IndexedTensor> for PartialTensor {
    /// Transfer an explicitly indexed tensor to a consumer supporting partial
    /// interfaces, preserving its payload, layout and established proof handles.
    fn from(value: IndexedTensor) -> Self {
        value.into_partial_engine()
    }
}

impl PartialTensor {
    /// Require explicit public indices without changing payload or assigning labels.
    ///
    /// Unchecked sources undergo the existing interface/scope validation once.
    /// Already established facts survive the representation-only conversion.
    pub fn try_into_indexed(self) -> Result<IndexedTensor> {
        self.ensure_validated()?;
        let Self {
            structure,
            expression,
            is_metric,
            is_composite,
            proofs,
        } = self;
        let structure = IndexedStructure::try_from(structure)
            .map_err(|error| TensorInferenceError::Invalid(error.to_string()))?;
        for (axis, slot) in structure.canonical().external_structure_iter().enumerate() {
            if crate::CookSettings::index_has_open_port(slot.aind)
                .map_err(|error| TensorInferenceError::Invalid(error.to_string()))?
            {
                return Err(TensorInferenceError::Invalid(
                    spenso::structure::StructureError::ExpectedExplicitSlot { axis }.to_string(),
                ));
            }
        }
        Ok(IndexedTensor {
            structure,
            expression,
            is_metric,
            is_composite,
            proofs,
        })
    }
}

impl IndexedTensor {
    /// Admit source syntax and require every free tensor port to be explicitly indexed.
    pub fn infer_indexed(expression: Atom) -> Result<Self> {
        PartialTensor::infer(expression)?.try_into_indexed()
    }

    /// Check a declared indexed boundary, including the logical shape of zero.
    /// Encoded unresolved ports are rejected even though `AbstractIndex` can store them.
    pub fn from_indexed_parts(expression: Atom, structure: IndexedStructure) -> Result<Self> {
        let is_metric = expression
            .as_fun_view()
            .is_some_and(|fun| fun.get_symbol() == spenso::network::library::symbolic::ETS.metric);
        let is_composite = !matches!(expression.as_view(), AtomView::Fun(_));
        // Observe these declared parts directly. Rebuilding through logical-slot
        // sorting could change canonical tie order for symbolic dimensions and
        // discard the supplied layout. The checked conversion below validates
        // callbacks, scopes and multiplicities before publishing this carrier.
        PartialTensor {
            structure: structure.into_partial(),
            expression,
            is_metric,
            is_composite,
            proofs: Default::default(),
        }
        .try_into_indexed()
    }

    pub(crate) fn into_partial_engine(self) -> PartialTensor {
        let Self {
            structure,
            expression,
            is_metric,
            is_composite,
            proofs,
        } = self;
        PartialTensor {
            structure: structure.into_partial(),
            expression,
            is_metric,
            is_composite,
            proofs,
        }
    }

    fn borrowed_parts<'a>(&'a self, structure: &'a PartialStructure) -> BorrowedTensorParts<'a> {
        BorrowedTensorParts {
            expression: &self.expression,
            structure,
            is_metric: self.is_metric,
            is_composite: self.is_composite,
            proofs: &self.proofs,
        }
    }

    fn partial_engine(&self) -> PartialTensor {
        self.clone().into_partial_engine()
    }

    pub fn rank(&self) -> usize {
        self.structure.canonical().order()
    }

    pub fn is_scalar(&self) -> bool {
        self.rank() == 0
    }

    pub fn is_zero(&self) -> bool {
        self.expression.is_zero()
    }

    pub fn same_value(&self, other: &Self) -> bool {
        self.expression == other.expression
            && self.structure.logical_slots() == other.structure.logical_slots()
    }

    pub fn hash_value<H: Hasher>(&self, state: &mut H) {
        use std::hash::Hash;
        self.expression.hash(state);
        self.structure.logical_slots().hash(state);
    }

    pub fn contraction_complete(&self) -> bool {
        self.proofs.contracted
    }

    pub fn reduction_status(&self) -> ReductionStatus {
        self.proofs
            .reduction
            .as_ref()
            .map_or(ReductionStatus::Deferred, |record| record.status)
    }

    pub fn with_checked_expression(&self, expression: Atom) -> Result<Self> {
        let structure = self.structure.clone().into_partial();
        self.borrowed_parts(&structure)
            .with_checked_expression(expression)?
            .try_into_indexed()
    }

    pub fn with_algebra_result(&self, expression: Atom) -> Result<Self> {
        let structure = self.structure.clone().into_partial();
        self.borrowed_parts(&structure)
            .with_algebra_result(expression)?
            .try_into_indexed()
    }

    pub fn with_rewritten_expression(&self, expression: Atom) -> Result<Self> {
        let structure = self.structure.clone().into_partial();
        self.borrowed_parts(&structure)
            .with_rewritten_expression(expression)?
            .try_into_indexed()
    }

    pub fn with_transformed_expression(&self, expression: Atom) -> Result<(Self, bool)> {
        let structure = self.structure.clone().into_partial();
        let (value, contracted) = self
            .borrowed_parts(&structure)
            .with_transformed_expression(expression)?;
        Ok((value.try_into_indexed()?, contracted))
    }

    pub fn with_lorentz_dimension(&self, dimension: Dimension) -> Result<Self> {
        self.partial_engine()
            .with_lorentz_dimension(dimension)?
            .try_into_indexed()
    }

    pub fn map_representation_dimensions(
        &self,
        map: impl FnMut(LibraryRep, AtomView<'_>) -> Atom,
    ) -> Result<Self> {
        self.partial_engine()
            .map_representation_dimensions(map)?
            .try_into_indexed()
    }

    pub fn with_cooked_indices(&self, settings: &crate::CookSettings) -> Result<Self> {
        self.partial_engine()
            .with_cooked_indices(settings)?
            .try_into_indexed()
    }

    pub fn wrap_indices(&self, scope: Symbol, dummies_only: bool) -> Result<Self> {
        self.partial_engine()
            .wrap_indices(scope, dummies_only)?
            .try_into_indexed()
    }

    pub fn rename_symbols(&self, renaming: &HashMap<Symbol, Symbol>) -> Result<Self> {
        self.partial_engine()
            .rename_symbols(renaming)?
            .try_into_indexed()
    }

    pub fn expanded(&self, var: Option<AtomView<'_>>, via_poly: bool) -> Result<Self> {
        self.partial_engine()
            .expanded(var, via_poly)?
            .try_into_indexed()
    }

    pub fn try_add(&self, right: &Self) -> Result<Self> {
        self.partial_engine()
            .try_add(&right.partial_engine())?
            .try_into_indexed()
    }

    pub fn try_sub(&self, right: &Self) -> Result<Self> {
        self.partial_engine()
            .try_sub(&right.partial_engine())?
            .try_into_indexed()
    }

    pub fn pow(&self, exponent: &Self) -> Result<Self> {
        self.partial_engine()
            .pow(&exponent.partial_engine())?
            .try_into_indexed()
    }

    pub fn divide(&self, denominator: &Self) -> Result<Self> {
        self.partial_engine()
            .divide(&denominator.partial_engine())?
            .try_into_indexed()
    }

    /// Consume both operands through the shared composition engine and strict
    /// indexed result validation. Clone an operand explicitly when retaining it.
    pub fn multiply(self, right: Self) -> std::result::Result<Self, TensorCompositionError> {
        let left = self.into_partial_engine();
        let right = right.into_partial_engine();
        left.multiply(&right)?
            .try_into_indexed()
            .map_err(composition_error)
    }

    pub fn reindex_interface_ports(
        &self,
        replacements: &HashMap<usize, AbstractIndex>,
    ) -> std::result::Result<Self, TensorCompositionError> {
        self.partial_engine()
            .reindex_interface_ports(replacements)?
            .try_into_indexed()
            .map_err(composition_error)
    }

    pub fn simplify_algebra(&self, settings: &AlgebraSettings) -> Result<Self> {
        self.partial_engine()
            .into_plan_reduction(ReductionRequest::Algebra(settings))?
            .try_into_indexed()
    }

    pub fn contract(&self, settings: ContractSettings<'_>) -> Result<Self> {
        self.partial_engine()
            .into_plan_reduction(ReductionRequest::Contract(settings))?
            .try_into_indexed()
    }

    pub fn canonize(&self) -> std::result::Result<Self, CanonicalizationError> {
        Ok(self.partial_engine().canonize()?.try_into_indexed()?)
    }

    pub fn to_dots(&self) -> Result<Self> {
        self.partial_engine().to_dots()?.try_into_indexed()
    }

    pub fn expand_dots(&self) -> Result<Self> {
        self.partial_engine().expand_dots()?.try_into_indexed()
    }

    pub fn undo_chain(&self) -> Result<Self> {
        self.partial_engine().undo_chain()?.try_into_indexed()
    }

    pub fn undo_trace(&self) -> Result<Self> {
        self.partial_engine().undo_trace()?.try_into_indexed()
    }

    pub fn collect<const N: usize>(&self, filter: TensorCollectFilter<N>) -> Result<Self> {
        self.partial_engine().collect(filter)?.try_into_indexed()
    }

    pub fn coefficients_equal<const N: usize>(
        &self,
        right: &Self,
        filter: TensorCollectFilter<N>,
    ) -> Result<ConditionResult> {
        self.partial_engine()
            .coefficients_equal(&right.partial_engine(), filter)
    }

    pub fn coefficients_are_zero<const N: usize>(
        &self,
        filter: TensorCollectFilter<N>,
    ) -> Result<ConditionResult> {
        self.partial_engine().coefficients_are_zero(filter)
    }

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
                    .map(|slot| slot.aind),
                &state,
            );
        }
        state
    }

    pub fn dirac_adjoint(&self, preserve_indices: bool) -> std::result::Result<Self, AdjointError> {
        Ok(self
            .partial_engine()
            .dirac_adjoint(preserve_indices)?
            .try_into_indexed()?)
    }

    /// Consume the indexed carrier at the existing finite-component boundary.
    /// The original logical layout and ranked-zero handling remain authoritative.
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn into_network<T, S, K, Aind, Lib, FunLib>(
        self,
        library: &Lib,
        functions: &FunLib,
        settings: &ParseSettings,
        index_decoding: Option<&crate::CookSettings>,
    ) -> std::result::Result<
        Canonicalized<Network<NetworkStore<T, Atom>, K, Symbol, Aind>>,
        TensorNetworkError<K, Symbol>,
    >
    where
        Aind: AbsInd + DummyAind + ParseableAind,
        K: Clone + std::fmt::Debug + std::fmt::Display,
        S: TensorStructure<Slot = LibrarySlot<Aind>>
            + ScalarStructure
            + Clone
            + StructureFromAtom
            + From<OrderedStructure<LibraryRep, Aind>>,
        T: HasStructure<Structure = S>
            + TensorStructure<Slot = LibrarySlot<Aind>>
            + Clone
            + From<SparseTensor<Atom, S>>
            + for<'src> TensorFromExpression<'src, S, Atom, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Atom, Key = Symbol>,
    {
        self.into_partial_engine()
            .into_network(library, functions, settings, index_decoding)
    }
}

fn composition_error(error: TensorInferenceError) -> TensorCompositionError {
    TensorCompositionError::InvalidResultInterface(error.to_string())
}

impl std::hash::Hash for IndexedTensor {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.hash_value(state);
    }
}

#[cfg(feature = "bincode")]
impl Encode for IndexedTensor {
    fn encode<E: Encoder>(&self, encoder: &mut E) -> StdResult<(), EncodeError> {
        let structure = self.structure.clone().into_partial();
        self.borrowed_parts(&structure).encode(encoder)
    }
}

#[cfg(feature = "bincode")]
impl<C: HasStateMap> Decode<C> for IndexedTensor {
    fn decode<D: Decoder<Context = C>>(decoder: &mut D) -> StdResult<Self, DecodeError> {
        PartialTensor::decode(decoder)?
            .try_into_indexed()
            .map_err(|error| DecodeError::OtherString(error.to_string()))
    }
}

#[cfg(test)]
mod tests {

    #[cfg(feature = "bincode")]
    use super::super::composition::ExplicitIndexOccurrences;
    #[cfg(feature = "bincode")]
    use super::super::inference::tests::{
        INFERENCE_CALLS, OCCURRENCE_FUNCTION_VISITS, SCOPE_VALIDATIONS,
    };
    #[cfg(feature = "bincode")]
    use super::super::persistence::TENSOR_ARCHIVE_VERSION;
    use super::*;
    use crate::{CookSettings, test_support};
    #[cfg(feature = "bincode")]
    use bincode::config;
    use spenso::iterators::IteratableTensor;
    use spenso::network::library::panicing::ErroringLibrary;
    use spenso::network::library::symbolic::{ExplicitKey, TensorLibrary};
    use spenso::network::library::{DummyKey, DummyLibrary};
    use spenso::network::parsing::ShadowedStructure;
    use spenso::network::{ExecutionResult, Sequential, SmallestDegree};
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    use spenso::structure::partial::{PartialIndex, PartialStructureExt};
    use spenso::structure::representation::{Euclidean, Minkowski, RepName};
    use spenso::structure::slot::IsAbstractSlot;
    use spenso::tensors::parametric::ParamTensor;
    use spenso::{index_symbol, tensor_symbol};
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    #[cfg(feature = "bincode")]
    use symbolica::state::StateMap;
    use symbolica::{function, symbol};

    #[test]
    fn indexed_admission_keeps_layout_and_rejects_all_unresolved_ports() {
        crate::test_support::test_initialize();
        let rep = Minkowski {}.new_rep(4).to_lib();
        let slots = [
            rep.slot(AbstractIndex::Normal(19)),
            rep.slot(AbstractIndex::Normal(7)),
        ];
        let structure: IndexedStructure = slots.into_iter().collect();
        let zero = IndexedTensor::from_indexed_parts(Atom::Zero, structure.clone()).unwrap();
        assert_eq!(zero.structure.logical_slots(), slots);
        assert_eq!(zero.structure.layout(), structure.layout());
        assert_eq!(zero.rank(), 2);
        assert!(!zero.is_scalar());
        assert!(zero.is_zero());
        let restored = zero
            .clone()
            .into_partial_engine()
            .try_into_indexed()
            .unwrap();
        assert_eq!(restored, zero);
        assert!(restored.same_reduction_state(&zero));

        let head = spenso::tensor_symbol!("indexed_strict_open_vector");
        assert!(IndexedTensor::infer_indexed(function!(head, rep.to_symbolic([]))).is_err());
        let open = AbstractIndex::Open { owner: 31, axis: 0 };
        for index in [open, open.scoped(symbol!("indexed_strict_scope"))] {
            let interface: IndexedStructure = [rep.slot(index)].into_iter().collect();
            assert!(IndexedTensor::from_indexed_parts(Atom::Zero, interface).is_err());
            let partial = PartialTensor::new(
                Atom::Zero,
                PartialStructure::from_logical_slots([rep.slot(PartialIndex::Explicit(index))]),
            );
            assert!(partial.try_into_indexed().is_err());
        }
        let malformed = PartialTensor::new(Atom::one(), structure.into_partial());
        assert!(malformed.try_into_indexed().is_err());

        // A symbolic dimension can compare equal to a concrete dimension during
        // canonical sorting. A metadata projection must not sort those ties anew.
        let symbolic = Minkowski {}.new_rep(symbol!("indexed_layout_D")).to_lib();
        let tied: IndexedStructure = [
            symbolic.slot(AbstractIndex::Normal(23)),
            rep.slot(AbstractIndex::Normal(23)),
        ]
        .into_iter()
        .collect();
        let partial = PartialTensor {
            expression: Atom::Zero,
            structure: tied.clone().into_partial(),
            is_metric: false,
            is_composite: true,
            proofs: Default::default(),
        };
        let indexed = partial.try_into_indexed().unwrap();
        assert_eq!(indexed.structure, tied);
        assert_eq!(indexed.structure.layout(), tied.layout());
    }

    #[test]
    fn indexed_metadata_projection_keeps_completed_proofs_and_graph_identity() {
        crate::test_support::test_initialize();
        let source = PartialTensor::infer(symbolica::parse!("indexed_proof_x + 1"))
            .unwrap()
            .simplify_algebra(&AlgebraSettings::default())
            .unwrap();
        assert_eq!(source.reduction_status(), ReductionStatus::Complete);
        let graph = source.shallow_graph().unwrap();
        let indexed = source.clone().try_into_indexed().unwrap();
        let roundtrip: PartialTensor = indexed.clone().into();
        assert!(source.same_reduction_state(&roundtrip));
        assert!(std::sync::Arc::ptr_eq(
            &graph,
            &roundtrip.shallow_graph().unwrap()
        ));
        let identity = indexed
            .map_representation_dimensions(|_, dim| dim.to_owned())
            .unwrap();
        assert!(indexed.same_reduction_state(&identity));
        let fixed = identity
            .simplify_algebra(&AlgebraSettings::default())
            .unwrap();
        assert!(identity.same_reduction_state(&fixed));
        assert_eq!(fixed.reduction_status(), ReductionStatus::Complete);
    }

    #[test]
    fn owned_completed_reduction_retains_payload_allocation_and_cached_graph() {
        crate::test_support::test_initialize();
        let settings = AlgebraSettings::default();
        let source = PartialTensor::infer(symbolica::parse!(
            "(owned_reduction_x+owned_reduction_y)^17"
        ))
        .unwrap()
        .simplify_algebra(&settings)
        .unwrap();
        assert_eq!(source.reduction_status(), ReductionStatus::Complete);
        let graph = source.shallow_graph().unwrap();
        let expected = source.clone();
        // The shared borrowed entry still leaves the caller's allocation intact;
        // its owned counterpart can reuse that allocation at completed return.
        let borrowed = source
            .plan_reduction(ReductionRequest::Algebra(&settings))
            .unwrap();
        assert!(borrowed.same_value(&source));
        assert!(borrowed.same_reduction_state(&source));
        let storage = source.expression().as_view().get_data().as_ptr();
        let result = source
            .into_plan_reduction(ReductionRequest::Algebra(&settings))
            .unwrap();
        assert_eq!(result.expression().as_view().get_data().as_ptr(), storage);
        assert!(result.same_value(&expected));
        assert!(result.same_reduction_state(&expected));
        assert!(std::sync::Arc::ptr_eq(
            &graph,
            &result.shallow_graph().unwrap()
        ));
    }

    #[test]
    fn indexed_checked_publication_borrows_source_proof_cells_and_consumes_result() {
        crate::test_support::test_initialize();
        let slot = Euclidean {}
            .new_rep(2)
            .to_lib()
            .slot::<AbstractIndex, _>(211);
        let expression = function!(
            spenso::tensor_symbol!("borrowed_publication_vector"),
            slot.to_atom()
        );
        let mut source = IndexedTensor::infer_indexed(expression).unwrap();
        source.proofs.validated = Default::default();
        source.proofs.algebra = Default::default();
        let result = source
            .with_algebra_result(Atom::num(2) * source.expression())
            .unwrap();
        assert_eq!(source.proofs.algebra.get(), Some(&true));
        assert_eq!(source.proofs.validated.get(), Some(&true));
        assert_eq!(result.structure.logical_slots(), vec![slot]);
        let structure = source.structure.clone().into_partial();
        let borrowed = source.borrowed_parts(&structure);
        assert!(std::ptr::eq(borrowed.expression, source.expression()));
        assert!(std::ptr::eq(borrowed.proofs, &source.proofs));

        let replacement = Atom::num(3) * source.expression();
        let allocation = replacement.as_view().get_data().as_ptr();
        let checked = source.with_checked_expression(replacement).unwrap();
        assert_eq!(
            checked.expression().as_view().get_data().as_ptr(),
            allocation
        );
        assert_eq!(checked.structure.logical_slots(), vec![slot]);
        assert!(source.with_checked_expression(Atom::one()).is_err());
        let zero = source.with_rewritten_expression(Atom::Zero).unwrap();
        assert_eq!(zero.structure.logical_slots(), vec![slot]);
    }

    #[test]
    fn indexed_dummy_reservation_reuses_borrowed_observations_and_public_ports() {
        use super::super::composition::DUMMY_RESERVATION_WALKS;
        crate::test_support::test_initialize();
        let rep = Minkowski {}.new_rep(4).to_lib();
        let public = AbstractIndex::Normal(223);
        let hidden = AbstractIndex::Normal(227);
        let metadata = function!(
            symbol!("indexed_reservation_metadata"; Scalar),
            rep.slot::<AbstractIndex, _>(hidden).to_atom()
        );
        let source = PartialTensor::infer(
            metadata
                * function!(
                    spenso::tensor_symbol!("indexed_reservation_vector"),
                    rep.slot::<AbstractIndex, _>(public).to_atom()
                ),
        )
        .unwrap();
        let observations = std::sync::Arc::clone(source.reduction_observations());
        assert!(observations.reserved_indices.contains(&hidden));
        let indexed = source.clone().try_into_indexed().unwrap();
        DUMMY_RESERVATION_WALKS.with(|count| count.set(0));
        let partial_names = PartialTensor::reserved_dummies([&source]);
        let indexed_names = IndexedTensor::reserved_dummies([&indexed]);
        assert_eq!(DUMMY_RESERVATION_WALKS.with(|count| count.get()), 0);
        assert!(std::sync::Arc::ptr_eq(
            indexed.proofs.observations.get().unwrap(),
            &observations
        ));
        for index in [partial_names.fresh_index(), indexed_names.fresh_index()] {
            assert_ne!(index, public);
            assert_ne!(index, hidden);
        }

        // Public labels are serialized in logical order. A named label's
        // callback can observe that order even though reservation uses a set.
        let calls = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let observed_calls = std::sync::Arc::clone(&calls);
        let name = symbol!(
            "indexed_reservation_order",
            tag = &spenso::network::tags::SPENSO_TAG.index,
            norm = move |node, _| {
                if let AtomView::Fun(function) = node {
                    observed_calls
                        .lock()
                        .unwrap()
                        .push(usize::try_from(function.iter().nth(1).unwrap()).unwrap());
                }
            }
        );
        let logical = [
            rep.slot(AbstractIndex::Named(name.into(), 0, 19)),
            rep.slot(AbstractIndex::Named(name.into(), 0, 7)),
        ];
        let indexed =
            IndexedTensor::from_indexed_parts(Atom::Zero, logical.into_iter().collect()).unwrap();
        assert_ne!(indexed.structure.canonical().external_structure(), logical);
        let partial: PartialTensor = indexed.clone().into();
        calls.lock().unwrap().clear();
        let _ = PartialTensor::reserved_dummies([&partial]);
        assert_eq!(*calls.lock().unwrap(), [19, 7]);
        calls.lock().unwrap().clear();
        let _ = IndexedTensor::reserved_dummies([&indexed]);
        assert_eq!(*calls.lock().unwrap(), [19, 7]);
    }

    #[test]
    fn indexed_admission_rejects_cooked_open_ports_without_rewriting_explicit_labels() {
        test_support::test_initialize();
        let cooking = CookSettings::reversible().preserve_tags();
        let rep = Minkowski {}.new_rep(4).to_lib();
        let scope = symbol!("indexed_cooked_open_scope");
        let open = AbstractIndex::Open {
            owner: 173,
            axis: 1,
        };
        for source in [open, open.scoped(scope)] {
            let encoded = cooking.try_cook(source.to_atom().as_view()).unwrap();
            let encoded = rep.to_symbolic([encoded]);
            let slot = LibrarySlot::<AbstractIndex>::try_from(encoded.as_view()).unwrap();
            assert!(matches!(slot.aind, AbstractIndex::Symbol(_)));
            let nested = function!(
                AbstractIndex::scope_symbol(),
                Atom::var(scope),
                slot.aind.to_atom()
            );
            let nested = cooking
                .cook_function(nested.as_fun_view().unwrap())
                .unwrap();
            let nested = AbstractIndex::try_from(nested.as_view()).unwrap();
            for index in [slot.aind, slot.aind.scoped(scope), nested] {
                let slot = rep.slot(index);
                let interface: IndexedStructure = [slot].into_iter().collect();
                assert!(IndexedTensor::from_indexed_parts(Atom::Zero, interface.clone()).is_err());
                let nonzero = PartialTensor::infer(function!(
                    tensor_symbol!("indexed_encoded_open_vector"),
                    slot.to_atom()
                ))
                .unwrap();
                assert!(nonzero.try_into_indexed().is_err());
            }
        }

        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let arbitrary = symbol!(
            "indexed_opaque_label",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let named = index_symbol!("indexed_named_label");
        for label in [
            function!(arbitrary, 17),
            function!(named, 7, 9),
            function!(AIND_SYMBOLS.cind, 0),
            function!(arbitrary, open.to_atom()),
        ] {
            let encoded = cooking.try_cook(label.as_view()).unwrap();
            let encoded = rep.to_symbolic([encoded]);
            let slot = LibrarySlot::<AbstractIndex>::try_from(encoded.as_view()).unwrap();
            assert!(matches!(slot.aind, AbstractIndex::Symbol(_)));
            let source = PartialTensor::from_declared_parts(
                Atom::Zero,
                PartialStructure::from_logical_slots([rep.slot(PartialIndex::Explicit(slot.aind))]),
            )
            .unwrap()
            .simplify_algebra(&AlgebraSettings::default())
            .unwrap();
            calls.store(0, Ordering::Relaxed);
            let indexed = source.clone().try_into_indexed().unwrap();
            let restored: PartialTensor = indexed.into();
            assert_eq!(restored.expression(), source.expression());
            assert_eq!(restored.structure(), source.structure());
            assert!(source.same_reduction_state(&restored));
            assert_eq!(calls.load(Ordering::Relaxed), 0);
        }
    }

    #[test]
    fn indexed_conversion_does_not_replay_callbacks_and_checks_changed_results() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = spenso::tensor_symbol!(
            "indexed_conversion_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let slot = Euclidean {}
            .new_rep(2)
            .to_lib()
            .slot(AbstractIndex::Normal(47));
        let source = PartialTensor::infer(function!(head, slot.to_atom())).unwrap();
        calls.store(0, Ordering::Relaxed);
        let indexed = source.try_into_indexed().unwrap();
        let same = indexed
            .with_checked_expression(indexed.expression.clone())
            .unwrap();
        assert_eq!(calls.load(Ordering::Relaxed), 0);
        assert!(same.same_value(&indexed));
        let open = function!(
            spenso::tensor_symbol!("indexed_result_open"),
            slot.rep().to_symbolic([])
        );
        assert!(indexed.with_transformed_expression(open).is_err());
        assert!(indexed.with_checked_expression(Atom::one()).is_err());
        let zero = indexed.with_algebra_result(Atom::Zero).unwrap();
        assert_eq!(zero.structure.logical_slots(), vec![slot]);
        let invalid = spenso::tensor_symbol!(
            "indexed_bad_callback",
            norm = |_, output| {
                **output = Atom::one();
            }
        );
        assert!(
            IndexedTensor::from_indexed_parts(
                function!(invalid, slot.to_atom()),
                [slot].into_iter().collect()
            )
            .is_err()
        );
        assert!(
            IndexedTensor::from_indexed_parts(
                Atom::var(spenso::network::tags::SPENSO_TAG.chain_in),
                std::iter::empty::<LibrarySlot<AbstractIndex>>().collect()
            )
            .is_err()
        );
    }

    #[test]
    fn indexed_algebra_and_dimension_mapping_keep_the_declared_boundary() {
        crate::test_support::test_initialize();
        let dimension = symbol!("indexed_dimension_D");
        let rep = Minkowski {}.new_rep(dimension).to_lib();
        let slots = [
            rep.slot(AbstractIndex::Normal(71)),
            rep.slot(AbstractIndex::Normal(69)),
        ];
        let expression = function!(
            spenso::network::library::symbolic::ETS.metric,
            slots[0].to_atom(),
            slots[1].to_atom()
        );
        let source =
            IndexedTensor::from_indexed_parts(expression, slots.into_iter().collect()).unwrap();
        let reduced = source
            .simplify_algebra(&AlgebraSettings::default())
            .unwrap();
        assert_eq!(reduced.structure.logical_slots(), slots);
        for value in [
            reduced.clone(),
            reduced.with_algebra_result(Atom::Zero).unwrap(),
        ] {
            let mapped = value
                .map_representation_dimensions(|_, _| Atom::num(4))
                .unwrap();
            let actual = mapped.structure.logical_slots();
            assert_eq!(
                actual.iter().map(|slot| slot.aind).collect::<Vec<_>>(),
                vec![slots[0].aind, slots[1].aind]
            );
            assert!(
                actual
                    .iter()
                    .all(|slot| slot.dim() == Dimension::Concrete(4))
            );
            assert_eq!(mapped.is_zero(), value.is_zero());
        }
    }

    #[test]
    fn indexed_component_lowering_retains_nonsymmetric_logical_axes_and_ranked_zero() {
        crate::test_support::test_initialize();
        type Structure = ShadowedStructure<AbstractIndex>;
        type Tensor = ParamTensor<Structure>;
        type Lib = TensorLibrary<ParamTensor<ExplicitKey<AbstractIndex>>, AbstractIndex>;
        type Net =
            Network<NetworkStore<Tensor, Atom>, ExplicitKey<AbstractIndex>, Symbol, AbstractIndex>;
        let reps = [
            Minkowski {}.new_rep(4).to_lib(),
            Euclidean {}.new_rep(2).to_lib(),
        ];
        let slots = [
            reps[0].slot(AbstractIndex::Normal(89)),
            reps[1].slot(AbstractIndex::Normal(83)),
        ];
        let head = spenso::tensor_symbol!("indexed_component_matrix");
        let source = IndexedTensor::from_indexed_parts(
            function!(head, slots[0].to_atom(), slots[1].to_atom()),
            slots.into_iter().collect(),
        )
        .unwrap();
        let mut library = Lib::new();
        library
            .insert_explicit_sparse(
                ExplicitKey::from_iter(reps, head, None),
                (0..4).flat_map(|row| {
                    (0..2).map(move |col| (vec![row, col], Atom::num(10 * row + col + 1)))
                }),
                Atom::Zero,
            )
            .unwrap();
        for (value, zero) in [
            (source.clone(), false),
            (source.with_algebra_result(Atom::Zero).unwrap(), true),
        ] {
            let lowered: Canonicalized<Net> = value
                .into_network(
                    &library,
                    &ErroringLibrary::<Symbol>::new(),
                    &ParseSettings::default(),
                    None,
                )
                .unwrap();
            let (mut network, layout) = lowered.into_parts();
            network
                .execute::<Sequential, SmallestDegree, _, _, _>(
                    &library,
                    &ErroringLibrary::<Symbol>::new(),
                )
                .unwrap();
            let ExecutionResult::Val(tensor) = network.result_tensor(&library).unwrap() else {
                panic!("matrix must remain a tensor, including zero");
            };
            assert_eq!(tensor.order(), 2);
            assert_eq!(
                layout.canonical_to_logical(&tensor.external_structure()),
                slots
            );
            let mut count = 0;
            for (indices, coefficient) in tensor.iter_expanded() {
                let logical = layout.canonical_to_logical(&indices);
                let expected = if zero {
                    Atom::Zero
                } else {
                    Atom::num(10 * logical[0] + logical[1] + 1)
                };
                assert_eq!(coefficient.to_owned(), expected);
                count += 1;
            }
            assert!(zero || count == 8);
        }
        // Scalar lowering shares the same consuming engine entry point.
        let scalar = IndexedTensor::infer_indexed(Atom::one()).unwrap();
        let _: Canonicalized<Network<NetworkStore<Tensor, Atom>, DummyKey, Symbol, AbstractIndex>> =
            scalar
                .into_network(
                    &DummyLibrary::<Tensor>::new(),
                    &ErroringLibrary::<Symbol>::new(),
                    &ParseSettings::default(),
                    None,
                )
                .unwrap();
    }

    #[cfg(feature = "bincode")]
    #[test]
    fn indexed_archive_keeps_zero_layout_and_rejects_unresolved_or_invalid_parts() {
        test_support::test_initialize();
        let rep = Euclidean {}.new_rep(2).to_lib();
        let slots = [
            rep.slot(AbstractIndex::Normal(103)),
            rep.slot(AbstractIndex::Normal(101)),
        ];
        let source =
            IndexedTensor::from_indexed_parts(Atom::Zero, slots.into_iter().collect()).unwrap();
        let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
        let (restored, used): (IndexedTensor, _) = bincode::decode_from_slice_with_context(
            &bytes,
            config::standard(),
            StateMap::default(),
        )
        .unwrap();
        assert_eq!(used, bytes.len());
        assert_eq!(restored, source);
        assert_eq!(restored.structure.logical_slots(), slots);
        assert!(restored.proofs.shallow_graph.get().is_none());
        let cooked_open = CookSettings::reversible()
            .preserve_tags()
            .try_cook(
                AbstractIndex::Open { owner: 47, axis: 1 }
                    .to_atom()
                    .as_view(),
            )
            .unwrap();
        let cooked_open = rep.to_symbolic([cooked_open]);
        let cooked_open = LibrarySlot::<AbstractIndex>::try_from(cooked_open.as_view()).unwrap();
        assert!(matches!(cooked_open.aind, AbstractIndex::Symbol(_)));
        for (expression, slots) in [
            (Atom::one(), slots.to_vec()),
            (
                Atom::Zero,
                vec![rep.slot(AbstractIndex::Open { owner: 43, axis: 0 })],
            ),
            (Atom::Zero, vec![cooked_open]),
        ] {
            let structure = slots
                .iter()
                .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind)))
                .collect::<Vec<_>>();
            let metadata = if expression.is_zero() {
                ExplicitIndexOccurrences::from_interface(&PartialStructure::from_logical_slots(
                    structure.iter().copied(),
                ))
            } else {
                ExplicitIndexOccurrences::default()
            };
            let bytes = bincode::encode_to_vec(
                (
                    TENSOR_ARCHIVE_VERSION,
                    expression,
                    structure,
                    false,
                    true,
                    None::<bool>,
                    metadata,
                ),
                config::standard(),
            )
            .unwrap();
            assert!(
                bincode::decode_from_slice_with_context::<_, IndexedTensor, _>(
                    &bytes,
                    config::standard(),
                    StateMap::default()
                )
                .is_err()
            );
        }
    }

    #[cfg(feature = "bincode")]
    #[test]
    fn indexed_archive_retains_cancelled_dummy_history_without_raw_readmission() {
        test_support::test_initialize();
        let slot = Euclidean {}
            .new_rep(3)
            .to_lib()
            .slot(AbstractIndex::Normal(97121));
        let x = Atom::var(symbol!("indexed_archive_cancel_x"));
        let polynomial = (&x + 1).pow(2) - &x * &x - Atom::num(2) * &x - 1;
        let source = PartialTensor::infer(
            function!(
                tensor_symbol!("indexed_archive_cancel_pair"),
                slot.to_atom(),
                slot.to_atom()
            ) * polynomial,
        )
        .unwrap()
        .expanded(None, false)
        .unwrap()
        .try_into_indexed()
        .unwrap();
        assert!(source.is_zero());
        let fresh =
            IndexedTensor::from_indexed_parts(source.expression.clone(), source.structure.clone())
                .unwrap();
        assert!(fresh.proofs.occurrences.get().unwrap().is_empty());
        let extra = IndexedTensor::infer_indexed(function!(
            tensor_symbol!("indexed_archive_cancel_extra"),
            slot.to_atom()
        ))
        .unwrap();
        INFERENCE_CALLS.with(|count| count.set(0));
        OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
        SCOPE_VALIDATIONS.with(|count| count.set(0));
        let bytes = bincode::encode_to_vec(&source, config::standard()).unwrap();
        let (restored, used): (IndexedTensor, _) = bincode::decode_from_slice_with_context(
            &bytes,
            config::standard(),
            StateMap::default(),
        )
        .unwrap();
        assert_eq!(used, bytes.len());
        assert!(restored.same_value(&source));
        let cloned = restored.clone();
        assert!(Arc::ptr_eq(
            restored.proofs.occurrences.get().unwrap(),
            cloned.proofs.occurrences.get().unwrap(),
        ));
        let partial: PartialTensor = cloned.into();
        let converted = partial.try_into_indexed().unwrap();
        for value in [source, restored, converted] {
            assert_eq!(value.proofs.occurrences.get().unwrap().count(slot), 2);
            assert!(matches!(
                value.multiply(extra.clone()),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            ));
        }
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
        assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
    }

    #[cfg(feature = "bincode")]
    #[test]
    fn indexed_archive_remaps_payload_dimension_and_public_index_symbols() {
        use base64::{Engine, engine::general_purpose::STANDARD};
        use std::io::Cursor;
        use symbolica::state::State;
        const CHILD: &str = "IDENSO_INDEXED_ARCHIVE_CHILD";
        let tensor = || {
            let dimension = symbol!("idenso::indexed_archive::dimension");
            let index = symbol!("idenso::indexed_archive::index");
            IndexedTensor::infer_indexed(function!(
                spenso::tensor_symbol!("idenso::indexed_archive::vector"),
                Minkowski {}
                    .new_rep(dimension)
                    .to_symbolic([Atom::var(index)])
            ))
            .unwrap()
        };
        if std::env::var_os(CHILD).is_some() {
            let source = tensor();
            let bytes = bincode::encode_to_vec(source, bincode::config::standard()).unwrap();
            let mut state = Vec::new();
            State::export(&mut state).unwrap();
            println!(
                "INDEXED_INDEX={}",
                symbol!("idenso::indexed_archive::index").get_id()
            );
            println!("INDEXED_STATE={}", STANDARD.encode(state));
            println!("INDEXED_BYTES={}", STANDARD.encode(bytes));
            return;
        }
        let output = std::process::Command::new(std::env::current_exe().unwrap()).args([
            "--exact", "tensor::indexed::tests::indexed_archive_remaps_payload_dimension_and_public_index_symbols", "--nocapture"
        ]).env(CHILD, "1").output().unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let stdout = String::from_utf8(output.stdout).unwrap();
        let field = |prefix| {
            stdout
                .lines()
                .find_map(|line| line.strip_prefix(prefix))
                .unwrap()
        };
        let old_index: u32 = field("INDEXED_INDEX=").parse().unwrap();
        for position in 0..=old_index {
            let _ = symbol!(&format!("idenso::indexed_archive::padding_{position}"));
        }
        let expected = tensor();
        assert_ne!(
            symbol!("idenso::indexed_archive::index").get_id(),
            old_index
        );
        let map = State::import(
            &mut Cursor::new(STANDARD.decode(field("INDEXED_STATE=")).unwrap()),
            None,
        )
        .unwrap();
        let bytes = STANDARD.decode(field("INDEXED_BYTES=")).unwrap();
        let (restored, used): (IndexedTensor, _) =
            bincode::decode_from_slice_with_context(&bytes, bincode::config::standard(), map)
                .unwrap();
        assert_eq!(used, bytes.len());
        assert_eq!(restored, expected);
        assert_eq!(
            restored.structure.logical_slots(),
            expected.structure.logical_slots()
        );
    }
}
