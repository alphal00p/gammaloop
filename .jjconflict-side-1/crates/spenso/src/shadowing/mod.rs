use std::env;

use crate::{
    network::{library::function_lib::INBUILTS, library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        CanonicalLayout, HasName, HasStructure, TensorShell, TensorStructure, ToSymbolic,
        abstract_index::AIND_SYMBOLS,
        concrete_index::FlatIndex,
        slot::{IsAbstractSlot, ParseableAind},
    },
    symbolica_init::in_symbolica_initializer,
    tensors::{
        data::{DataTensor, DenseTensor},
        parametric::{ExpandedCoefficent, MixedTensor, ParamTensor, TensorCoefficient},
        // symbolic::SymbolicTensor,
    },
};
use ::symbolica_utils::{IntoArgs, IntoSymbol};
use eyre::Result;
use symbolica::{atom::Atom, evaluate::FunctionMap, initialize};

mod atom_conversion;
mod collect;
mod macros;
mod projectors;
pub mod static_symbols;
mod trace;

pub mod symbolica_utils;

pub use atom_conversion::IntoAtom;
pub use collect::{
    COLLECT, Collectable, TensorCollectExt, TensorCollectFilter, TermLeaf, TermTape,
};
pub(crate) use projectors::expand_chain_like_projector;
pub use projectors::{ANTISYM, CYCLIC, ProjectorExpander, SYM, antisym, cyclic, sym};
pub use trace::{trace, trace_factor_views, trace_parts, trace_sym};

initialize!(|| {
    in_symbolica_initializer(|| {
        let _ = INBUILTS.force_in_initializer().conj;
        let _ = *ANTISYM.force_in_initializer();
        let _ = *CYCLIC.force_in_initializer();
        let _ = *SYM.force_in_initializer();
        let _ = *COLLECT.force_in_initializer();
        let _ = ETS.force_in_initializer().delta;
        let _ = ETS.force_in_initializer().metric;
        let _ = SPENSO_TAG.force_in_initializer().bracket;
        let _ = AIND_SYMBOLS.force_in_initializer().scope;
    });
});

#[cfg(test)]
mod tests;

/// Trait that enables shadowing of a tensor
///
/// This creates a dense tensor of atoms, where the atoms are the expanded indices of the tensor, with the global name as the name of the labels.
pub trait Shadowable:
    HasStructure<
        Structure: TensorStructure + HasName<Name: IntoSymbol, Args: IntoArgs> + Clone + Sized,
    > + Sized
where
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    // type Const;
    fn expanded_shadow(&self) -> Result<DenseTensor<Atom, Self::Structure>> {
        self.shadow(Self::Structure::expanded_coef)
    }
    fn expanded_shadow_logical(
        &self,
        layout: &CanonicalLayout,
    ) -> Result<DenseTensor<Atom, Self::Structure>> {
        self.shadow(|structure, id| {
            let index = structure.co_expanded_index(id).unwrap();
            let index = layout.canonical_to_logical(&index).into();
            ExpandedCoefficent {
                name: structure.name().map(|name| name.ref_into_symbol()),
                index,
                args: structure.args(),
            }
        })
    }

    fn flat_shadow(&self) -> Result<DenseTensor<Atom, Self::Structure>> {
        self.shadow(Self::Structure::flat_coef)
    }

    fn shadow<T>(
        &self,
        index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> T,
    ) -> Result<DenseTensor<Atom, Self::Structure>>
    where
        T: TensorCoefficient,
    {
        self.structure().clone().to_dense_labeled(index_to_atom)
    }
}

pub trait Concretize<T> {
    /// Materialize the target representation, returning errors from finite-component
    /// construction. Symbolic targets may retain symbolic dimensions.
    fn concretize(self) -> Result<T>;
    fn concretize_logical(self, layout: &CanonicalLayout) -> Result<T>;
}

fn sparse_shadow_tensors() -> bool {
    env::var_os("SPENSO_SPARSE_SHADOW_TENSORS").is_some()
}

impl<S: Shadowable> Concretize<DenseTensor<Atom, S::Structure>> for S
where
    <<S::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn concretize(self) -> Result<DenseTensor<Atom, S::Structure>> {
        // self.flat_s
        // todo!()
        self.expanded_shadow()
    }

    fn concretize_logical(
        self,
        layout: &CanonicalLayout,
    ) -> Result<DenseTensor<Atom, S::Structure>> {
        self.expanded_shadow_logical(layout)
    }
}

impl<S: Shadowable> Concretize<DataTensor<Atom, S::Structure>> for S
where
    <<S::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn concretize(self) -> Result<DataTensor<Atom, S::Structure>> {
        // self.flat_s
        // todo!()
        let dense = <S as Concretize<DenseTensor<Atom, S::Structure>>>::concretize(self)?;
        Ok(if sparse_shadow_tensors() {
            dense.to_sparse().into()
        } else {
            dense.into()
        })
    }

    fn concretize_logical(
        self,
        layout: &CanonicalLayout,
    ) -> Result<DataTensor<Atom, S::Structure>> {
        let dense =
            <S as Concretize<DenseTensor<Atom, S::Structure>>>::concretize_logical(self, layout)?;
        Ok(if sparse_shadow_tensors() {
            dense.to_sparse().into()
        } else {
            dense.into()
        })
    }
}

impl<S: Shadowable> Concretize<ParamTensor<S::Structure>> for S
where
    <<S::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn concretize(self) -> Result<ParamTensor<S::Structure>> {
        // self.flat_s
        // todo!()
        <S as Concretize<DataTensor<Atom, S::Structure>>>::concretize(self).map(ParamTensor::param)
    }

    fn concretize_logical(self, layout: &CanonicalLayout) -> Result<ParamTensor<S::Structure>> {
        <S as Concretize<DataTensor<Atom, S::Structure>>>::concretize_logical(self, layout)
            .map(ParamTensor::param)
    }
}

impl<T: Clone, S: Shadowable> Concretize<MixedTensor<T, S::Structure>> for S
where
    <<S::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn concretize(self) -> Result<MixedTensor<T, S::Structure>> {
        // self.flat_s
        // todo!()
        <S as Concretize<DataTensor<Atom, S::Structure>>>::concretize(self).map(MixedTensor::param)
    }

    fn concretize_logical(self, layout: &CanonicalLayout) -> Result<MixedTensor<T, S::Structure>> {
        <S as Concretize<DataTensor<Atom, S::Structure>>>::concretize_logical(self, layout)
            .map(MixedTensor::param)
    }
}

pub trait ShadowMapping: Shadowable
where
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn expanded_shadow_with_map(
        &self,
        fn_map: &mut FunctionMap,
    ) -> Result<ParamTensor<Self::Structure>> {
        self.shadow_with_map(fn_map, Self::Structure::expanded_coef)
    }

    fn shadow_with_map<T, F>(
        &self,
        fn_map: &mut FunctionMap,
        index_to_atom: F,
    ) -> Result<ParamTensor<Self::Structure>>
    where
        T: TensorCoefficient,
        F: Fn(&Self::Structure, FlatIndex) -> T + Clone,
    {
        // Some(ParamTensor::Param(self.shadow(index_to_atom)?.into()))
        self.append_map(fn_map, index_to_atom.clone());
        self.shadow(index_to_atom)
            .map(|x| ParamTensor::param(x.into()))
    }

    fn append_map<T>(
        &self,
        fn_map: &mut FunctionMap,
        index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> T,
    ) where
        T: TensorCoefficient;

    fn flat_append_map(&self, fn_map: &mut FunctionMap) {
        self.append_map(fn_map, Self::Structure::flat_coef)
    }

    fn expanded_append_map(&self, fn_map: &mut FunctionMap) {
        self.append_map(fn_map, Self::Structure::expanded_coef)
    }

    fn flat_shadow_with_map(
        &self,
        fn_map: &mut FunctionMap,
    ) -> Result<ParamTensor<Self::Structure>> {
        self.shadow_with_map(fn_map, Self::Structure::flat_coef)
    }
}

impl<S: TensorStructure + HasName + Clone> Shadowable for TensorShell<S>
where
    S::Name: IntoSymbol + Clone,
    S::Args: IntoArgs,
    <<S as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
}

impl<S: TensorStructure + HasName + Clone> ShadowMapping for TensorShell<S>
where
    S::Name: IntoSymbol + Clone,
    S::Args: IntoArgs,
    <<S as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn append_map<T>(
        &self,
        _fn_map: &mut FunctionMap,
        _index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> T,
    ) where
        T: TensorCoefficient,
    {
    }
}

#[cfg(test)]
pub mod test {
    use std::sync::RwLock;

    use once_cell::sync::Lazy;
    use symbolica::atom::Symbol;

    use crate::network::library::{symbolic::ETS, symbolic::ExplicitKey, symbolic::TensorLibrary};

    #[allow(clippy::type_complexity)]
    pub static EXPLICIT_TENSOR_MAP: Lazy<
        RwLock<TensorLibrary<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>>,
    > = Lazy::new(|| {
        let mut lib = TensorLibrary::new();
        lib.update_ids();
        RwLock::new(lib)
    });

    use crate::structure::abstract_index::AbstractIndex;
    use crate::structure::{Canonicalized, OrderedStructure};
    use crate::{
        contraction::Contract,
        structure::{
            HasStructure, IndexlessNamedStructure, TensorStructure,
            representation::{LibraryRep, RepName},
        },
        tensors::parametric::MixedTensor,
    };

    #[test]
    fn test_identity() {
        // EXPLICIT_TENSOR_MAP.write().unwrap().update_ids();
        let mut tensor_library: TensorLibrary<
            MixedTensor<f64, ExplicitKey<AbstractIndex>>,
            AbstractIndex,
        > = TensorLibrary::new();
        tensor_library.update_ids();

        for rep in LibraryRep::all_representations() {
            let structure = [rep.new_rep(4), rep.new_rep(4).dual()];

            let idstructure: Canonicalized<IndexlessNamedStructure<Symbol, (), LibraryRep>> =
                IndexlessNamedStructure::from_iter(structure, ETS.metric, None);

            let idkey = ExplicitKey::from_structure(&idstructure).unwrap();

            let id = tensor_library
                .get(&idkey)
                .unwrap()
                .into_owned()
                .into_canonical();

            let trace_structure: OrderedStructure =
                Canonicalized::from_iter([rep.new_rep(4).slot(3), rep.new_rep(4).dual().slot(4)])
                    .into_canonical();
            let id1 = id.map_structure(|_| trace_structure.clone());
            let id2 = id1
                .clone()
                .map_structure(|_| trace_structure.clone().dual());

            // println!("{}", rep);
            assert_eq!(
                4.,
                id1.contract(&id2)
                    .unwrap()
                    .scalar()
                    .unwrap()
                    .try_into_concrete()
                    .unwrap()
                    .try_into_real()
                    .unwrap(),
                "trace of 4-dim identity should be 4 for rep {}",
                rep
            );
        }
    }

    // #[test]
    // fn pslash() {
    //     let expr = parse!("p(1,mink(4,mu))*γ(mink(4,mu),bis(4,i),bis(4,j))").unwrap();
    //     let mut network: TensorNetwork<MixedTensor<f64, ShadowedStructure>, Atom> =
    //         TensorNetwork::<MixedTensor<f64, ShadowedStructure>, Atom>::try_from_view(
    //             expr.as_view(),
    //             &EXPLICIT_TENSOR_MAP.read().unwrap(),
    //         )
    //         .unwrap();
    //     network.contract().unwrap();
    //     let result = network
    //         .result()
    //         .unwrap()
    //         .0
    //         .try_into_parametric()
    //         .unwrap()
    //         .tensor
    //         .map_data(|a| a.to_string())
    //         .map_structure(VecStructure::from);

    //     insta::assert_ron_snapshot!(result);
    // }
}
