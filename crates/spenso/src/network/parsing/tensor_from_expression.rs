use std::fmt::{Debug, Display};
use std::ops::AddAssign;

use symbolica::atom::{Atom, AtomOrView, AtomView, Symbol};

use super::{LeafInterfaceCache, ParseSettings, ShorthandParsing, StructureFromAtom};
use crate::algebra::ScalarMul;
use crate::network::library::{FunctionLibrary, Library, LibraryTensor};
use crate::network::store::NetworkStore;
use crate::network::{
    ContractionStrategy, ExecuteOp, ExecutionResult, FastTensorSum, MinIntermediateCost, Network,
    Ref, Sequential, TensorNetworkError,
};
use crate::shadowing::Concretize;
use crate::structure::permuted::ApplyPendingIndexPermutation;
use crate::structure::slot::{IsAbstractSlot, SlotMatcher};
use crate::structure::{
    Canonicalized, HasStructure, ScalarStructure, ScalarTensor, StructureError, TensorShell,
    TensorStructure,
};
use crate::tensors::complex::RealOrComplexTensor;
use crate::tensors::data::{DataTensor, DenseTensor, SparseTensor};
use crate::tensors::parametric::{EvalTensor, ParamOrConcrete, ParamTensor};

/// Builds a tensor leaf from an opaque expression once its structure is known.
///
/// Opaque shorthand parsing separates topology from realization: structure
/// inference decides the exposed slots, then this trait decides how the target
/// tensor type represents the original expression. Symbolic tensors can keep
/// the expression directly, while concretized tensor types may parse and
/// execute an expanded sub-network here with the MinIntermediateCost policy.
/// The returned layout describes the realized tensor: eager contraction can
/// remove internal ports from the originally inferred interface.
pub trait TensorFromExpression<'src, S, Sc, K, FK, Aind, Lib, FunLib>: Sized
where
    S: TensorStructure,
    Self: HasStructure,
    Lib: TensorLibraryFor<S, Self, Key = K>,
    FunLib: FunctionLibrary<Self, Sc, Key = FK>,
{
    /// Metadata retained by the source admission owner. Component tensors do
    /// not need an extra context; symbolic tensors reuse their regional facts.
    type Admission: Clone;

    /// Read a previously admitted boundary. Targets without an admission owner
    /// reject this path instead of falling back to structural inference.
    fn admitted_structure(
        _expression: AtomView<'_>,
        _admission: &Self::Admission,
        _next_index: impl FnMut() -> Aind,
    ) -> Result<Canonicalized<S>, StructureError> {
        Err(StructureError::ParsingError(
            "tensor target has no admitted structure provider".into(),
        ))
    }

    /// Share operation-local recognition state across admitted parser leaves.
    fn admitted_structure_with_matcher(
        expression: AtomView<'_>,
        admission: &Self::Admission,
        _matcher: &mut SlotMatcher,
        next_index: impl FnMut() -> Aind,
    ) -> Result<Canonicalized<S>, StructureError> {
        Self::admitted_structure(expression, admission, next_index)
    }

    #[allow(clippy::result_large_err)]
    fn tensor_from_expression(
        expression: AtomOrView<'src>,
        structure: Canonicalized<S>,
        tensor_library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        admission: Option<&Self::Admission>,
        leaf_interfaces: Option<&LeafInterfaceCache>,
    ) -> Result<Canonicalized<Self>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display;

    /// Realize an ordinary tensor leaf after the library has declined it.
    /// Unlike an opaque composite, this leaf must not recursively parse itself.
    #[allow(clippy::result_large_err)]
    fn tensor_from_leaf(
        expression: AtomOrView<'src>,
        structure: Canonicalized<S>,
        admission: Option<&Self::Admission>,
    ) -> Result<Canonicalized<Self>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display;
}

/// A tensor library whose stored tensor can be indexed into `T`.
pub trait TensorLibraryFor<S, T>:
    Library<S, Value = Canonicalized<<Self as TensorLibraryFor<S, T>>::LibraryTensor>>
where
    S: TensorStructure,
    T: HasStructure,
{
    type LibraryTensor: LibraryTensor<WithIndices = T> + Clone;
}

impl<S, T, L, LT> TensorLibraryFor<S, T> for L
where
    S: TensorStructure,
    T: HasStructure,
    L: Library<S, Value = Canonicalized<LT>>,
    LT: LibraryTensor<WithIndices = T> + Clone,
{
    type LibraryTensor = LT;
}

/// Marker for tensor types whose opaque expression can be realized by
/// expanding and executing a nested network.
pub trait ExpandedTensorFromExpression {}

impl<S> ExpandedTensorFromExpression for ParamTensor<S> {}
impl<C, S> ExpandedTensorFromExpression for ParamOrConcrete<C, S> {}
impl<T, S> ExpandedTensorFromExpression for EvalTensor<T, S> {}
impl<T, S: TensorStructure> ExpandedTensorFromExpression for RealOrComplexTensor<T, S> {}
impl<T, S> ExpandedTensorFromExpression for DataTensor<T, S> {}
impl<T, S> ExpandedTensorFromExpression for DenseTensor<T, S> {}
impl<T, S> ExpandedTensorFromExpression for SparseTensor<T, S> {}

impl<'src, S, Sc, T, K, Aind, Lib, FunLib>
    TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib> for T
where
    S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
    TensorShell<S>: Concretize<T>,
    S::Slot: IsAbstractSlot<Aind = Aind>,
    T::Slot: IsAbstractSlot<Aind = Aind>,
    T: ExpandedTensorFromExpression
        + HasStructure<Structure = S>
        + TensorStructure
        + Clone
        + ScalarTensor
        + ApplyPendingIndexPermutation<Output = T>
        + Ref
        + FastTensorSum
        + ScalarMul<Sc, Output = T>
        + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
    Sc: for<'r> TryFrom<AtomView<'r>> + TryFrom<Atom> + Clone + Into<T::Scalar>,
    TensorNetworkError<K, Symbol>: From<<Sc as TryFrom<Atom>>::Error>,
    for<'r> TensorNetworkError<K, Symbol>: From<<Sc as TryFrom<AtomView<'r>>>::Error>,
    K: Display + Debug + Clone,
    Aind: crate::structure::slot::AbsInd
        + crate::structure::slot::DummyAind
        + crate::structure::slot::ParseableAind,
    Lib: TensorLibraryFor<S, T, Key = K> + Sync,
    FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    NetworkStore<T, Sc>: ExecuteOp<FunLib, Lib, K, Symbol, Aind, Tensor = T, Scalar = Sc>,
    MinIntermediateCost: ContractionStrategy<NetworkStore<T, Sc>, Lib, K, Symbol, Aind>,
{
    type Admission = ();

    fn tensor_from_leaf(
        _expression: AtomOrView<'src>,
        structure: Canonicalized<S>,
        _admission: Option<&Self::Admission>,
    ) -> Result<Canonicalized<Self>, TensorNetworkError<K, Symbol>> {
        let (canonical, layout) = structure.into_parts();
        let tensor = canonical.to_shell().concretize_logical(&layout)?;
        Ok(Canonicalized::from_parts(tensor, layout))
    }

    fn tensor_from_expression(
        expression: AtomOrView<'src>,
        _structure: Canonicalized<S>,
        tensor_library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        _admission: Option<&Self::Admission>,
        leaf_interfaces: Option<&LeafInterfaceCache>,
    ) -> Result<Canonicalized<Self>, TensorNetworkError<K, Symbol>>
    where
        K: Display,
        Symbol: Display,
    {
        let mut expanded_settings = settings.clone();
        expanded_settings.shorthand_parsing = ShorthandParsing::expand_all();
        expanded_settings.depth_limit = None;

        let mut network = if let Some(leaf_interfaces) = leaf_interfaces {
            Network::<NetworkStore<T, Sc>, K, Symbol, Aind>::try_from_admitted_view_with_function_library::<S, Lib, FunLib>(
                expression.as_view(), tensor_library, function_library, &expanded_settings,
                leaf_interfaces.clone(), None,
            )?
        } else {
            Network::<NetworkStore<T, Sc>, K, Symbol, Aind>::try_from_view_with_function_library::<
                S,
                Lib,
                FunLib,
            >(
                expression.as_view(),
                tensor_library,
                function_library,
                &expanded_settings,
            )?
        };
        network.execute::<Sequential, MinIntermediateCost, Lib::LibraryTensor, Lib, FunLib>(
            tensor_library,
            function_library,
        )?;

        let tensor = match network.result_tensor(tensor_library)? {
            ExecutionResult::Val(tensor) => tensor.into_owned(),
            ExecutionResult::One => {
                let one: Sc = Atom::num(1).as_view().try_into()?;
                T::new_scalar(one.into())
            }
            ExecutionResult::Zero => {
                let zero: Sc = Atom::num(0).as_view().try_into()?;
                T::new_scalar(zero.into())
            }
        };
        // Execution has already applied logical permutations and contracted
        // internal ports. Its stored coordinate order is the interface to sew.
        Ok(Canonicalized::identity(tensor))
    }
}

#[cfg(test)]
mod tests {
    use symbolica::{
        atom::{Atom, AtomCore, Symbol},
        function, symbol,
    };

    use super::*;
    use crate::{
        network::library::{DummyKey, DummyLibrary, panicing::ErroringLibrary},
        structure::{
            HasStructure, NamedStructure,
            abstract_index::AbstractIndex,
            representation::{LibraryRep, Minkowski, RepName},
            slot::IsAbstractSlot,
        },
        tensor_symbol,
        tensors::parametric::ParamTensor,
    };

    #[test]
    fn param_tensor_from_expression_uses_inferred_structure() {
        type Structure = NamedStructure<Symbol, Vec<Atom>, LibraryRep, AbstractIndex>;

        let rep = Minkowski {}.new_rep(4);
        let slot = rep.slot(AbstractIndex::from(1));
        let structure: Canonicalized<Structure> =
            NamedStructure::from_iter([slot], symbol!("f"), None::<Vec<Atom>>);
        let expression = function!(
            tensor_symbol!(param_tensor_inferred_structure_opaque),
            slot.to_atom()
        );
        type TensorLib = DummyLibrary<ParamTensor<Structure>, DummyKey>;
        type FunLib = ErroringLibrary<Symbol>;
        let tensor = <ParamTensor<Structure> as TensorFromExpression<
            '_,
            Structure,
            Atom,
            DummyKey,
            Symbol,
            AbstractIndex,
            TensorLib,
            FunLib,
        >>::tensor_from_expression(
            expression.as_view().into(),
            structure,
            &TensorLib::new(),
            &FunLib::new(),
            &ParseSettings::default(),
            None,
            None,
        )
        .unwrap();

        assert_eq!(tensor.canonical().structure().order(), 1);
    }

    #[test]
    fn realized_opaque_chain_keeps_its_contracted_interface() {
        type Structure = NamedStructure<Symbol, Vec<Atom>, LibraryRep, AbstractIndex>;
        type Tensor = ParamTensor<Structure>;
        type Net = Network<NetworkStore<Tensor, Atom>, DummyKey, Symbol, AbstractIndex>;
        let rep = Minkowski {}.new_rep(2);
        let left = rep.slot::<AbstractIndex, _>(AbstractIndex::from(1));
        let right = rep.slot::<AbstractIndex, _>(AbstractIndex::from(2));
        let spectator = rep.slot::<AbstractIndex, _>(AbstractIndex::from(3));
        let tags = crate::network::tags::SPENSO_TAG;
        let factor = function!(
            tensor_symbol!(realized_opaque_matrix),
            tags.chain_in,
            tags.chain_out,
            spectator.to_atom()
        );
        // A weighted word may finish as a lazy ScaledTensor rather than a
        // LocalTensor. Eager realization must retain that scale too.
        let word = crate::chain!(left.to_atom(), right.to_atom(); [Atom::num(7) * &factor, factor]);
        let source = word
            * function!(tensor_symbol!(realized_opaque_left), left.to_atom())
            * function!(tensor_symbol!(realized_opaque_right), right.to_atom());
        let library = DummyLibrary::<Tensor>::new();
        let functions = ErroringLibrary::<Symbol>::new();
        let mut results = Vec::new();
        for shorthand_parsing in [
            ShorthandParsing::expand_all(),
            ShorthandParsing::expand_schoonschip_only(),
            ShorthandParsing::Opaque,
        ] {
            let mut network = Net::try_from_view::<Structure, _>(
                source.as_view(),
                &library,
                &ParseSettings {
                    shorthand_parsing,
                    ..Default::default()
                },
            )
            .unwrap();
            assert!(network.state.is_scalar());
            assert!(network.graph.dangling_indices().is_empty());
            network
                .execute::<Sequential, MinIntermediateCost, _, _, _>(&library, &functions)
                .unwrap();
            let scalar: Atom = network.result_scalar().unwrap().into();
            results.push(scalar);
        }
        assert_ne!(results[0], Atom::Zero);
        // Compare the finite two-dimensional component polynomial independently
        // of contraction schedule. Distinct vectors detect endpoint transposes.
        assert_eq!(results[0].expand(), results[1].expand());
        assert_eq!(results[0].expand(), results[2].expand());
    }
}
