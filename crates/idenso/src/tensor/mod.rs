use std::ops::AddAssign;

use linnet::half_edge::subgraph::subset::SubSet;
use spenso::{
    algebra::ScalarMul,
    bracket,
    contraction::{Contract, ContractionError, Trace},
    iterators::IteratableTensor,
    network::{
        ExecutionResult, Network, Ref, Sequential, SmallestDegree, TensorNetworkError,
        library::{
            DummyKey, DummyLibrary, FunctionLibrary, FunctionLibraryError,
            function_lib::Wrap,
            symbolic::{ETS, ExplicitKey, TensorLibrary},
        },
        parsing::{
            ParseSettings, ShadowedStructure, StrictTensorFilter, StructureFromAtom,
            StructureInferenceMode, TensorFromExpression, TensorLibraryFor,
        },
        store::NetworkStore,
        tags::SPENSO_TAG,
    },
    shadowing::{Concretize, symbolica_utils::SpensoPrintSettings},
    structure::{
        ApplyPendingIndexPermutation, CanonicalLayout, Canonicalized, HasName, HasStructure,
        MergeInfo, NamedStructure, OrderedStructure, PendingIndexPermutation, Reindexed,
        ScalarStructure, ScalarTensor, SlotIndex, StructureContract, TensorIdentity, TensorShell,
        TensorStructure, ToSymbolic,
        abstract_index::AIND_SYMBOLS,
        concrete_index::{ExpandedIndex, FlatIndex},
        representation::{LibraryRep, LibrarySlot},
        slot::{AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind},
    },
    tensors::parametric::MixedTensor,
};
use symbolica_utils::{AtomPrintExt, IntoArgs, IntoSymbol};

use delegate::delegate;
use spenso::structure::StructureError;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::dimension::Dimension;
use spenso::structure::representation::Representation;

use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    function,
};

use crate::NetworkToolingError;

pub mod aliases;
mod canonicalize;
pub mod composition;
pub mod inference;
mod replacement;
pub(crate) use canonicalize::remove_antisymmetric_zero_terms;
pub use replacement::TensorRule;

#[cfg(test)]
pub mod tests;

impl<S> FunctionLibrary<SymbolicTensor<S>, Atom> for Wrap {
    type Key = Symbol;
    fn apply(
        &self,
        key: &Self::Key,
        tensor: SymbolicTensor<S>,
    ) -> eyre::Result<SymbolicTensor<S>, FunctionLibraryError<Self::Key>> {
        Ok(SymbolicTensor {
            structure: tensor.structure,
            is_composite: true,
            is_metric: false,
            expression: function!(*key, tensor.expression),
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
        mode: StructureInferenceMode,
    ) -> Result<Canonicalized<Self>, StructureError> {
        let structure = OrderedStructure::<LibraryRep, Aind>::structure_from_atom(value, mode)?;
        Ok(structure.map_canonical(|structure| {
            let (is_composite, is_metric) = if let AtomView::Fun(f) = value {
                (false, f.get_symbol() == ETS.metric)
            } else {
                (true, false)
            };
            SymbolicTensor {
                structure,
                is_composite,
                is_metric,
                expression: value.to_owned(),
            }
        }))
    }
}

impl<Aind, K, Lib, FunLib>
    TensorFromExpression<
        SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
        Atom,
        K,
        Symbol,
        Aind,
        Lib,
        FunLib,
    > for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
where
    Aind: AbsInd + DummyAind + ParseableAind,
    K: std::fmt::Display,
    Lib: TensorLibraryFor<
            SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
            SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
            Key = K,
        >,
    FunLib: FunctionLibrary<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom, Key = Symbol>,
{
    fn tensor_from_expression(
        expression: AtomView<'_>,
        structure: Canonicalized<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>,
        _tensor_library: &Lib,
        _function_library: &FunLib,
        _settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        K: std::fmt::Display,
        Symbol: std::fmt::Display,
    {
        let mut tensor = structure.into_canonical();

        let (is_composite, is_metric) = if let AtomView::Fun(fun) = expression {
            (false, fun.get_symbol() == ETS.metric)
        } else {
            (true, false)
        };
        tensor.expression = expression.to_owned();
        tensor.is_composite = is_composite;
        tensor.is_metric = is_metric;
        Ok(tensor)
    }
}

/// A symbolic expression together with its tensor structure.
///
/// The default ordered, explicit-index structure supports network execution.
/// A [`spenso::structure::partial::PartialStructure`] also retains logical port
/// order and occurrence-local unresolved indices for positional composition.
/// The structure is independent of the expression, so a zero can retain its rank.
/// Checked construction and rewriting validate interfaces when callbacks can
/// change them; structure mappings retain the symbolic payload.
/// The payload defaults to `Atom`; borrowed views and Symbolica's `AliasedAtom`
/// use the same storage type. Aliased definitions retain their logical layouts in
/// [`aliases::AliasInterfaces`], separate from their single Symbolica registry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SymbolicTensor<S = OrderedStructure<LibraryRep, AbstractIndex>, E = Atom> {
    pub structure: S,
    pub is_metric: bool,
    pub is_composite: bool,
    pub expression: E,
}

impl<Aind: AbsInd + ParseableAind> SymbolicTensor<OrderedStructure<LibraryRep, Aind>> {
    /// Observe a callback-sensitive result without changing its planned port order.
    pub(crate) fn validate_rewritten_interface(
        &self,
        expression: &Atom,
    ) -> Result<(), inference::TensorInferenceError> {
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
                    inference::TensorInferenceError::Invalid(format!(
                        "cannot validate rewritten tensor index: {error}"
                    ))
                })?;
                if index.to_atom() != encoded {
                    return Err(inference::TensorInferenceError::Invalid(
                        "rewritten tensor index does not have a lossless shared representation"
                            .into(),
                    ));
                }
                Ok(slot.rep().slot(PartialIndex::Explicit(index)))
            })
            .collect::<Result<Vec<_>, inference::TensorInferenceError>>()?;
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
    fn new_scalar(scalar: Self::Scalar) -> Self {
        SymbolicTensor {
            structure: OrderedStructure::scalar_structure(),
            is_metric: false,
            is_composite: false,
            expression: scalar,
        }
    }
}

impl<Aind: AbsInd> ScalarStructure for SymbolicTensor<OrderedStructure<LibraryRep, Aind>> {
    fn scalar_structure() -> Self {
        SymbolicTensor {
            structure: OrderedStructure::scalar_structure(),
            is_metric: false,
            is_composite: false,
            expression: Atom::num(1),
        }
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

    fn apply_pending_index_permutation(
        mut self,
        pending: &PendingIndexPermutation,
    ) -> Self::Output {
        if pending.is_identity() {
            return self;
        }

        let target_slots = self.structure.external_structure();
        let mut dummy_structure = Vec::with_capacity(target_slots.len());
        let mut ids = Atom::one();
        for slot in pending.apply_slice_inverse(&target_slots) {
            let dummy = slot.to_dummy_ind();
            dummy_structure.push(dummy);
            ids *= Self::id(dummy.dual(), slot.to_lib()).expression;
        }

        let new_structure =
            Canonicalized::<OrderedStructure<LibraryRep, Aind>>::from_iter(dummy_structure)
                .into_canonical();

        for (o, n) in self
            .structure
            .external_structure_iter()
            .zip(new_structure.external_structure_iter())
        {
            self.expression = self.expression.replace(o.to_atom()).with(n.to_atom());
        }

        self.expression *= ids;
        self
    }
}

impl<Aind: AbsInd> TensorStructure for SymbolicTensor<OrderedStructure<LibraryRep, Aind>> {
    // type R = <T::Structure as TensorStructure>::R;
    type Indexed = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>;
    type Slot = LibrarySlot<Aind>;

    fn reindex_storage(self, indices: &[Aind]) -> Result<Reindexed<Self::Indexed>, StructureError> {
        let Self {
            structure,
            is_metric,
            is_composite,
            expression,
        } = self;
        Ok(structure
            .reindex_storage(indices)?
            .map_target(|structure| Self {
                structure,
                is_metric,
                is_composite,
                expression,
            }))
    }

    fn dual(self) -> Self {
        self.map_same_structure(|s| s.dual())
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

impl<S: TensorStructure> HasStructure for SymbolicTensor<S> {
    type Structure = S;
    type Scalar = Atom;
    type ScalarRef<'a>
        = &'a Atom
    where
        Self: 'a;
    type Store<O>
        = SymbolicTensor<O>
    where
        O: TensorStructure;

    fn map_structure<O: TensorStructure>(
        self,
        f: impl FnOnce(Self::Structure) -> O,
    ) -> Self::Store<O> {
        SymbolicTensor {
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
        &mut self.structure
    }

    fn map_same_structure(self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        SymbolicTensor {
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
        let expression = self.scalar()?;
        let mut indexed_product = false;
        if matches!(expression.as_view(), AtomView::Mul(_) | AtomView::Pow(_)) {
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
        }
        // Keep implicit index sums atomic when scalar arithmetic takes powers.
        Some(if indexed_product {
            bracket!(expression)
        } else {
            expression
        })
    }

    fn scalar_ref(&self) -> Option<&Self::Scalar> {
        if self.structure.is_scalar() {
            Some(&self.expression)
        } else {
            None
        }
    }
}

impl<Aind: AbsInd> IteratableTensor for SymbolicTensor<OrderedStructure<LibraryRep, Aind>> {
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

        Ok((
            Self {
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
        self.structure.trace_out();
    }

    fn trace(&mut self, i: usize, j: usize) {
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
        N: ToSymbolic + HasName + TensorStructure<Slot = LibrarySlot<Aind>>,
        N::Name: IntoSymbol + Clone,
        N::Args: IntoArgs,
    {
        let canonicalized = Canonicalized::from(structure.external_structure());
        let is_metric = structure.name()?.ref_into_symbol() == ETS.metric;
        Some(SymbolicTensor {
            expression: structure.to_symbolic(None)?,
            is_metric,
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
        let is_metric = canonical.name()?.ref_into_symbol() == ETS.metric;
        Some(SymbolicTensor {
            expression: structure.to_symbolic(None)?,
            is_composite: false,
            is_metric,
            structure: canonicalized.into_canonical(),
        })
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

        Network::<
            NetworkStore<MixedTensor<f64, ShadowedStructure<AbstractIndex>>, Atom>,
            ExplicitKey<AbstractIndex>,
            Symbol,
            AbstractIndex,
        >::try_from_view(self.expression.as_view(), library, &settings)
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

pub type SymbolicNet<Aind> = Network<
    NetworkStore<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom>,
    DummyKey,
    Symbol,
    Aind,
>;

// pub type HepTensor<Aind> = MixedTensor<f64, ShadowedStructure<Aind>>;

// pub type ParamNet<Aind> =
//     Network<NetworkStore<ParamTensor<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom>, DummyKey, Symbol, Aind>;
pub trait SymbolicNetExt<Aind: AbsInd + DummyAind + ParseableAind + 'static> {
    fn snapshot_dot(&self) -> String;
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

        SymbolicNet::<Aind>::try_from_view::<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, _>(
            *self, &lib, &settings,
        )
    }
}

pub mod contract;

impl<Aind: AbsInd + ParseableAind> Concretize<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
    for ShadowedStructure<Aind>
{
    fn concretize(self) -> eyre::Result<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>> {
        let is_metric = self.name() == Some(ETS.metric);
        Ok(SymbolicTensor {
            expression: self
                .to_symbolic(None)
                .ok_or_else(|| eyre::eyre!("symbolic tensor requires a name"))?,
            is_composite: false,
            is_metric,
            structure: self.structure,
        })
    }

    fn concretize_logical(
        self,
        layout: &CanonicalLayout,
    ) -> eyre::Result<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>> {
        let is_metric = self.name() == Some(ETS.metric);
        let logical_slots = layout.canonical_to_logical(&self.external_structure());
        let expression = FunctionBuilder::new(
            self.name()
                .ok_or_else(|| eyre::eyre!("symbolic tensor requires a name"))?,
        )
        .add_args(self.args().unwrap_or_default())
        .add_args(logical_slots.into_iter().map(|slot| slot.to_atom()))
        .finish();
        Ok(SymbolicTensor {
            expression,
            is_composite: false,
            is_metric,
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

impl<S> std::ops::Neg for SymbolicTensor<S> {
    type Output = SymbolicTensor<S>;

    fn neg(self) -> Self::Output {
        Self {
            expression: -self.expression,
            is_composite: true,
            is_metric: self.is_metric,
            structure: self.structure,
        }
    }
}

impl<Aind: AbsInd> AddAssign<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
    fn add_assign(&mut self, rhs: SymbolicTensor<OrderedStructure<LibraryRep, Aind>>) {
        debug_assert_eq!(self.structure, rhs.structure);
        self.expression += rhs.expression;
        self.is_composite = true;
        self.is_metric = false;
    }
}

impl<Aind: AbsInd> AddAssign<&SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
    for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
    fn add_assign(&mut self, rhs: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>) {
        debug_assert_eq!(self.structure, rhs.structure);
        self.expression += &rhs.expression;
        self.is_composite = true;
        self.is_metric = false;
    }
}

impl<S: Clone> ScalarMul<Atom> for SymbolicTensor<S> {
    type Output = SymbolicTensor<S>;

    fn scalar_mul(&self, rhs: &Atom) -> Option<Self::Output> {
        Some(Self {
            expression: rhs * &self.expression,
            is_composite: true,
            is_metric: false,
            structure: self.structure.clone(),
        })
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
