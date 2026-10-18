//! Fully indexed public boundaries backed by the existing tensor algebra.
//!
//! The engine projection only changes index metadata. It retains the payload,
//! logical layout, classification flags and immutable proof/cache handles. A
//! borrowed operation copies O(rank) metadata; it does not infer the expression
//! again or discard the engine's completed reduction state.

use std::{collections::HashMap, hash::Hasher};

use spenso::{
    network::{
        Network, TensorNetworkError,
        library::FunctionLibrary,
        parsing::{
            ParseSettings, ParseState, StructureFromAtom, TensorFromExpression, TensorLibraryFor,
        },
        store::NetworkStore,
    },
    shadowing::TensorCollectFilter,
    structure::{
        Canonicalized, HasStructure, IndexedStructure, OrderedStructure, ScalarStructure,
        TensorStructure,
        abstract_index::AbstractIndex,
        dimension::Dimension,
        partial::PartialStructure,
        representation::{LibraryRep, LibrarySlot},
        slot::{AbsInd, DummyAind, ParseableAind},
    },
    tensors::data::SparseTensor,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    id::ConditionResult,
};

use super::{
    AlgebraSettings, ContractSettings, SymbolicTensor, composition::TensorCompositionError,
    inference::TensorInferenceError, simplification::ReductionStatus,
};
use crate::{AdjointError, CanonicalizationError};

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
        self.partial_engine()
            .with_checked_expression(expression)?
            .try_into_indexed()
    }

    pub fn with_algebra_result(&self, expression: Atom) -> Result<Self> {
        self.partial_engine()
            .with_algebra_result(expression)?
            .try_into_indexed()
    }

    pub fn with_rewritten_expression(&self, expression: Atom) -> Result<Self> {
        self.partial_engine()
            .with_rewritten_expression(expression)?
            .try_into_indexed()
    }

    pub fn with_transformed_expression(&self, expression: Atom) -> Result<(Self, bool)> {
        let (value, contracted) = self
            .partial_engine()
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

    pub fn multiply(&self, right: &Self) -> std::result::Result<Self, TensorCompositionError> {
        self.partial_engine()
            .multiply(&right.partial_engine())?
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
            .simplify_algebra(settings)?
            .try_into_indexed()
    }

    pub fn contract(&self, settings: ContractSettings<'_>) -> Result<Self> {
        self.partial_engine().contract(settings)?.try_into_indexed()
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
        let values = values
            .into_iter()
            .map(Self::partial_engine)
            .collect::<Vec<_>>();
        PartialTensor::reserved_dummies(&values)
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
impl bincode::Encode for IndexedTensor {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> std::result::Result<(), bincode::error::EncodeError> {
        self.expression.encode(encoder)?;
        self.structure.logical_slots().encode(encoder)?;
        self.is_metric.encode(encoder)?;
        self.is_composite.encode(encoder)
    }
}

#[cfg(feature = "bincode")]
impl<C: symbolica::state::HasStateMap> bincode::Decode<C> for IndexedTensor {
    fn decode<D: bincode::de::Decoder<Context = C>>(
        decoder: &mut D,
    ) -> std::result::Result<Self, bincode::error::DecodeError> {
        use spenso::structure::representation::LibrarySlot;
        let expression = Atom::decode(decoder)?;
        let slots = Vec::<LibrarySlot<AbstractIndex>>::decode(decoder)?;
        let is_metric = bool::decode(decoder)?;
        let is_composite = bool::decode(decoder)?;
        let mut value = Self::from_indexed_parts(expression, slots.into_iter().collect())
            .map_err(|error| bincode::error::DecodeError::OtherString(error.to_string()))?;
        value.is_metric = is_metric;
        value.is_composite = is_composite;
        Ok(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::{
        iterators::IteratableTensor,
        network::{
            ExecutionResult, Sequential, SmallestDegree,
            library::{
                DummyKey, DummyLibrary,
                panicing::ErroringLibrary,
                symbolic::{ExplicitKey, TensorLibrary},
            },
            parsing::ShadowedStructure,
        },
        structure::{
            partial::{PartialIndex, PartialStructureExt},
            representation::{Euclidean, Minkowski, RepName},
            slot::IsAbstractSlot,
        },
        tensors::parametric::ParamTensor,
    };
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
    fn indexed_admission_rejects_cooked_open_ports_without_rewriting_explicit_labels() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        crate::test_support::test_initialize();
        let cooking = crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
        let rep = Minkowski {}.new_rep(4).to_lib();
        let scope = symbol!("indexed_cooked_open_scope");
        let open = AbstractIndex::Open {
            owner: 173,
            axis: 1,
        };
        for source in [open, open.scoped(scope)] {
            let encoded = cooking
                .try_cook_indices(rep.slot::<AbstractIndex, _>(source).to_atom().as_view())
                .unwrap();
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
                    spenso::tensor_symbol!("indexed_encoded_open_vector"),
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
        let named = spenso::index_symbol!("indexed_named_label");
        for label in [
            function!(arbitrary, 17),
            function!(named, 7, 9),
            function!(spenso::structure::abstract_index::AIND_SYMBOLS.cind, 0),
            function!(arbitrary, open.to_atom()),
        ] {
            let raw = function!(rep.rep.symbol(), 4, label);
            let encoded = cooking.try_cook_indices(raw.as_view()).unwrap();
            let slot = LibrarySlot::<AbstractIndex>::try_from(encoded.as_view()).unwrap();
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
        use symbolica::state::StateMap;
        crate::test_support::test_initialize();
        let rep = Euclidean {}.new_rep(2).to_lib();
        let slots = [
            rep.slot(AbstractIndex::Normal(103)),
            rep.slot(AbstractIndex::Normal(101)),
        ];
        let source =
            IndexedTensor::from_indexed_parts(Atom::Zero, slots.into_iter().collect()).unwrap();
        let bytes = bincode::encode_to_vec(&source, bincode::config::standard()).unwrap();
        let (restored, used): (IndexedTensor, _) = bincode::decode_from_slice_with_context(
            &bytes,
            bincode::config::standard(),
            StateMap::default(),
        )
        .unwrap();
        assert_eq!(used, bytes.len());
        assert_eq!(restored, source);
        assert_eq!(restored.structure.logical_slots(), slots);
        assert!(restored.proofs.shallow_graph.get().is_none());
        let cooked_open = crate::CookSettings::indices()
            .with_mode(crate::CookMode::ReversibleEncoding)
            .try_cook_indices(
                rep.slot::<AbstractIndex, _>(AbstractIndex::Open { owner: 47, axis: 1 })
                    .to_atom()
                    .as_view(),
            )
            .unwrap();
        let cooked_open = LibrarySlot::<AbstractIndex>::try_from(cooked_open.as_view()).unwrap();
        for (expression, slots) in [
            (Atom::one(), slots.to_vec()),
            (
                Atom::Zero,
                vec![rep.slot(AbstractIndex::Open { owner: 43, axis: 0 })],
            ),
            (Atom::Zero, vec![cooked_open]),
        ] {
            let bytes = bincode::encode_to_vec(
                (expression, slots, false, true),
                bincode::config::standard(),
            )
            .unwrap();
            assert!(
                bincode::decode_from_slice_with_context::<_, IndexedTensor, _>(
                    &bytes,
                    bincode::config::standard(),
                    StateMap::default()
                )
                .is_err()
            );
        }
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
