use std::fmt::{Debug, Display};

use symbolica::atom::{Atom, AtomView};

use super::{ParseSettings, ShorthandParsing, StructureFromAtom};
use crate::{
    network::{
        ContractionStrategy, ExecuteOp, ExecutionResult, Network, Sequential, SmallestDegree,
        TensorNetworkError, TensorOrScalarOrKey,
        library::{FunctionLibrary, Library, LibraryTensor},
        store::NetworkStore,
    },
    shadowing::Concretize,
    structure::{
        HasStructure, PermutedStructure, ScalarStructure, ScalarTensor, TensorShell,
        TensorStructure, permuted::PermuteTensor, slot::IsAbstractSlot,
    },
    tensors::{
        complex::RealOrComplexTensor,
        data::{DataTensor, DenseTensor, SparseTensor},
        parametric::{EvalTensor, ParamOrConcrete, ParamTensor},
    },
};

/// Builds a tensor leaf from an opaque expression once its structure is known.
///
/// Opaque shorthand parsing separates topology from realization: structure
/// inference decides the exposed slots, then this trait decides how the target
/// tensor type represents the original expression. Symbolic tensors can keep
/// the expression directly, while concretized tensor types may parse and
/// execute an expanded sub-network here.
pub trait TensorFromExpression<S, Sc, K, FK, Aind, Lib, FunLib>: Sized
where
    S: TensorStructure,
    Self: HasStructure,
    Lib: TensorLibraryFor<S, Self, Key = K>,
    FunLib: FunctionLibrary<Self, Sc, Key = FK>,
{
    #[allow(clippy::result_large_err)]
    fn tensor_from_expression(
        expression: AtomView<'_>,
        structure: PermutedStructure<S>,
        tensor_library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display;
}

/// A tensor library whose stored tensor can be indexed into `T`.
pub trait TensorLibraryFor<S, T>:
    Library<S, Value = PermutedStructure<<Self as TensorLibraryFor<S, T>>::LibraryTensor>>
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
    L: Library<S, Value = PermutedStructure<LT>>,
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

impl<S, Sc, T, K, Aind, Lib, FunLib>
    TensorFromExpression<S, Sc, K, symbolica::atom::Symbol, Aind, Lib, FunLib> for T
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
        + PermuteTensor<Permuted = T>,
    Sc: for<'r> TryFrom<AtomView<'r>> + Clone + Into<T::Scalar>,
    for<'r> TensorNetworkError<K, symbolica::atom::Symbol>:
        From<<Sc as TryFrom<AtomView<'r>>>::Error>,
    K: Display + Debug + Clone,
    Aind: crate::structure::slot::AbsInd
        + crate::structure::slot::DummyAind
        + crate::structure::slot::ParseableAind,
    Lib: TensorLibraryFor<S, T, Key = K> + Sync,
    FunLib: FunctionLibrary<T, Sc, Key = symbolica::atom::Symbol>,
    NetworkStore<T, Sc>: ExecuteOp<FunLib, Lib, K, symbolica::atom::Symbol, Aind>,
    SmallestDegree: ContractionStrategy<NetworkStore<T, Sc>, Lib, K, symbolica::atom::Symbol, Aind>,
{
    fn tensor_from_expression(
        expression: AtomView<'_>,
        _structure: PermutedStructure<S>,
        tensor_library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, symbolica::atom::Symbol>>
    where
        K: Display,
        symbolica::atom::Symbol: Display,
    {
        let mut expanded_settings = settings.clone();
        expanded_settings.shorthand_parsing = ShorthandParsing::expand_all();

        let mut network =
            Network::<NetworkStore<T, Sc>, K, symbolica::atom::Symbol, Aind>::try_from_view_with_function_library::<
                S,
                Lib,
                FunLib,
            >(
                expression,
                tensor_library,
                function_library,
                &expanded_settings,
            )?;
        network.execute::<Sequential, SmallestDegree, Lib::LibraryTensor, Lib, FunLib>(
            tensor_library,
            function_library,
        )?;

        match network.result()? {
            ExecutionResult::Val(TensorOrScalarOrKey::Tensor { tensor, .. }) => Ok(tensor.clone()),
            ExecutionResult::Val(TensorOrScalarOrKey::Scalar(scalar)) => {
                Ok(T::new_scalar(scalar.clone().into()))
            }
            ExecutionResult::Val(TensorOrScalarOrKey::Key { nodeid, .. }) => {
                Ok(network
                    .graph
                    .get_lib_data::<S, Lib::LibraryTensor, Lib>(tensor_library, nodeid)?)
            }
            ExecutionResult::One => {
                let one: Sc = Atom::num(1).as_view().try_into()?;
                Ok(T::new_scalar(one.into()))
            }
            ExecutionResult::Zero => {
                let zero: Sc = Atom::num(0).as_view().try_into()?;
                Ok(T::new_scalar(zero.into()))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use symbolica::{
        atom::{Atom, Symbol},
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
        let structure: PermutedStructure<Structure> =
            NamedStructure::from_iter([slot], symbol!("f"), None::<Vec<Atom>>);
        let expression = function!(tensor_symbol!(opaque), slot.to_atom());
        type TensorLib = DummyLibrary<ParamTensor<Structure>, DummyKey>;
        type FunLib = ErroringLibrary<Symbol>;
        let tensor = <ParamTensor<Structure> as TensorFromExpression<
            Structure,
            Atom,
            DummyKey,
            Symbol,
            AbstractIndex,
            TensorLib,
            FunLib,
        >>::tensor_from_expression(
            expression.as_view(),
            structure,
            &TensorLib::new(),
            &FunLib::new(),
            &ParseSettings::default(),
        )
        .unwrap();

        assert_eq!(tensor.structure().order(), 1);
    }
}
