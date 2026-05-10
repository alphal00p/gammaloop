use std::fmt::Display;

use symbolica::atom::AtomView;

use super::ParseSettings;
use crate::{
    network::{
        TensorNetworkError,
        library::{FunctionLibrary, Library},
    },
    shadowing::Concretize,
    structure::{PermutedStructure, TensorShell, TensorStructure},
};

/// Builds a tensor leaf from an opaque expression once its structure is known.
///
/// Opaque shorthand parsing separates topology from realization: structure
/// inference decides the exposed slots, then this trait decides how the target
/// tensor type represents the original expression. Symbolic tensors can keep
/// the expression directly, while concretized tensor types may parse and
/// execute an expanded sub-network here.
pub trait TensorFromExpression<S, Sc, K, FK, Aind>: Sized
where
    S: TensorStructure,
{
    #[allow(clippy::result_large_err)]
    fn tensor_from_expression<Lib, FunLib>(
        expression: AtomView<'_>,
        structure: PermutedStructure<S>,
        tensor_library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
        Lib: Library<S, Key = K>,
        FunLib: FunctionLibrary<Self, Sc, Key = FK>;
}

impl<S, Sc, T, K, FK, Aind> TensorFromExpression<S, Sc, K, FK, Aind> for T
where
    S: TensorStructure,
    TensorShell<S>: Concretize<T>,
{
    fn tensor_from_expression<Lib, FunLib>(
        _expression: AtomView<'_>,
        structure: PermutedStructure<S>,
        _tensor_library: &Lib,
        _function_library: &FunLib,
        _settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
        Lib: Library<S, Key = K>,
        FunLib: FunctionLibrary<Self, Sc, Key = FK>,
    {
        Ok(structure
            .structure
            .to_shell()
            .concretize(Some(structure.index_permutation.inverse())))
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
        tensors::parametric::ParamTensor,
    };

    #[test]
    fn param_tensor_from_expression_uses_inferred_structure() {
        type Structure = NamedStructure<Symbol, Vec<Atom>, LibraryRep, AbstractIndex>;

        let rep = Minkowski {}.new_rep(4);
        let slot = rep.slot(AbstractIndex::from(1));
        let structure: PermutedStructure<Structure> =
            NamedStructure::from_iter([slot], symbol!("f"), None::<Vec<Atom>>);
        let expression = function!(symbol!("opaque"), slot.to_atom());
        let tensor = <ParamTensor<Structure> as TensorFromExpression<
            Structure,
            Atom,
            DummyKey,
            Symbol,
            AbstractIndex,
        >>::tensor_from_expression(
            expression.as_view(),
            structure,
            &DummyLibrary::<ParamTensor<Structure>, DummyKey>::new(),
            &ErroringLibrary::<Symbol>::new(),
            &ParseSettings::default(),
        )
        .unwrap();

        assert_eq!(tensor.structure().order(), 1);
    }
}
