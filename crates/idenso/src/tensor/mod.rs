use std::ops::AddAssign;

use linnet::{
    half_edge::subgraph::{SubSetLike, subset::SubSet},
    permutation::Permutation,
};
use spenso::{
    algebra::ScalarMul,
    contraction::{Contract, ContractionError},
    network::{
        ExecutionResult, Network, Ref, Sequential, SmallestDegree, TensorNetworkError,
        TensorOrScalarOrKey,
        library::{
            DummyKey, DummyLibrary, FunctionLibrary, FunctionLibraryError, TensorLibraryData,
            function_lib::Wrap,
            symbolic::{ETS, ExplicitKey, TensorLibrary},
        },
        parsing::{Parse, ParseSettings, ShadowedStructure},
        store::NetworkStore,
    },
    shadowing::{
        Concretize,
        symbolica_utils::{AtomCoreExt, IntoArgs, IntoSymbol},
    },
    structure::{
        HasName, HasStructure, MergeInfo, NamedStructure, OrderedStructure, PermutedStructure,
        ScalarStructure, ScalarTensor, SlotIndex, StructureContract, TensorShell, TensorStructure,
        ToSymbolic,
        abstract_index::AIND_SYMBOLS,
        permuted::PermuteTensor,
        representation::{LibraryRep, LibrarySlot},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind},
    },
    tensors::parametric::MixedTensor,
};

use delegate::delegate;
use spenso::structure::StructureError;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::dimension::Dimension;
use spenso::structure::representation::Representation;

use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    id::Replacement,
    symbol,
};

#[cfg(test)]
pub mod tests;

impl<Aind: AbsInd> FunctionLibrary<SymbolicTensor<Aind>, Atom> for Wrap {
    type Key = Symbol;
    fn apply(
        &self,
        key: &Self::Key,
        tensor: SymbolicTensor<Aind>,
    ) -> eyre::Result<SymbolicTensor<Aind>, FunctionLibraryError<Self::Key>> {
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

impl<Aind: ParseableAind + AbsInd> Parse for SymbolicTensor<Aind> {
    fn parse(value: AtomView) -> Result<PermutedStructure<Self>, StructureError> {
        let structure = OrderedStructure::from_atomcore_unchecked(value)?;
        let (is_composite, is_metric) = if let AtomView::Fun(f) = value {
            (false, f.get_symbol() == ETS.metric)
        } else {
            (true, false)
        };
        Ok(PermutedStructure::identity(SymbolicTensor {
            structure,
            is_composite,
            is_metric,
            expression: value.to_owned(),
        }))
    }
}
/// A fully symbolic tensor, with no concrete values.
///
/// This tensor is used to represent the structure of a tensor, and is used to perform symbolic contraction.
/// Currently contraction is just a multiplication of the atoms, but in the future this will ensure that internal indices are independent accross the contraction.
///
/// Additionally, this can also be used as a tensor structure, that tracks the history, much like [`HistoryStructure`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SymbolicTensor<Aind: AbsInd = AbstractIndex> {
    pub structure: OrderedStructure<LibraryRep, Aind>,
    pub is_metric: bool,
    pub is_composite: bool,
    pub expression: symbolica::atom::Atom,
}

impl<Aind: AbsInd> Ref for SymbolicTensor<Aind> {
    type Ref<'a>
        = &'a SymbolicTensor<Aind>
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        self
    }
}

impl<Aind: AbsInd> ScalarTensor for SymbolicTensor<Aind> {
    fn new_scalar(scalar: Self::Scalar) -> Self {
        SymbolicTensor {
            structure: OrderedStructure::scalar_structure(),
            is_metric: false,
            is_composite: false,
            expression: scalar,
        }
    }
}

impl<Aind: AbsInd + DummyAind + ParseableAind> PermuteTensor for SymbolicTensor<Aind> {
    type Id = SymbolicTensor<Aind>;
    type IdSlot = LibrarySlot<Aind>;
    type Permuted = SymbolicTensor<Aind>;

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        Self::from_named(&NamedStructure::<Symbol, (), LibraryRep, Aind>::id(i, j)).unwrap()
    }

    fn permute_inds(mut self, permutation: &linnet::permutation::Permutation) -> Self::Permuted {
        let (new_structure, idstructures) = self.structure.clone().permute_inds(permutation);

        for (o, n) in self
            .structure
            .external_structure_iter()
            .zip(new_structure.external_structure_iter())
        {
            self.expression = self.expression.replace(o.to_atom()).with(n.to_atom());
        }

        let mut ids = Atom::one();
        for s in idstructures.iter() {
            let o = s.external_structure();
            ids *= Self::id(o[0], o[1]).expression;
        }
        self.expression *= ids;

        self
    }

    fn permute_reps(mut self, rep_perm: &linnet::permutation::Permutation) -> Self::Permuted {
        let (new_structure, idstructures) = self.structure.clone().permute_reps(rep_perm);

        for (o, n) in self
            .structure
            .external_structure_iter()
            .zip(new_structure.external_structure_iter())
        {
            self.expression = self.expression.replace(o.to_atom()).with(n.to_atom());
        }

        let mut ids = Atom::one();
        for s in idstructures.iter() {
            let o = s.external_structure();
            ids *= Self::id(o[0], o[1]).expression;
        }

        self.expression *= ids;

        self
    }
}

impl<Aind: AbsInd> TensorStructure for SymbolicTensor<Aind> {
    // type R = <T::Structure as TensorStructure>::R;
    type Indexed = SymbolicTensor<Aind>;
    type Slot = LibrarySlot<Aind>;

    fn reindex(self, indices: &[Aind]) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        let res = self.structure.reindex(indices)?;
        Ok(PermutedStructure {
            structure: Self {
                structure: res.structure,
                is_metric: self.is_metric,
                is_composite: self.is_composite,
                expression: self.expression,
            },
            rep_permutation: res.rep_permutation,
            index_permutation: res.index_permutation,
        })
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

impl<Aind: AbsInd> HasStructure for SymbolicTensor<Aind> {
    type Structure = OrderedStructure<LibraryRep, Aind>;
    type Scalar = Atom;
    type ScalarRef<'a>
        = &'a Atom
    where
        Self: 'a;
    type Store<S>
        = TensorShell<S>
    where
        S: TensorStructure;

    fn map_structure<O: TensorStructure>(
        self,
        f: impl FnOnce(Self::Structure) -> O,
    ) -> Self::Store<O> {
        TensorShell {
            structure: f(self.structure),
        }
    }

    fn map_structure_result<O: TensorStructure, Er>(
        self,
        f: impl FnOnce(Self::Structure) -> eyre::Result<O, Er>,
    ) -> std::result::Result<Self::Store<O>, Er> {
        Ok(TensorShell {
            structure: f(self.structure)?,
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
        if self.is_scalar() {
            Some(self.expression)
        } else {
            None
        }
    }

    fn scalar_ref(&self) -> Option<&Self::Scalar> {
        if self.is_scalar() {
            Some(&self.expression)
        } else {
            None
        }
    }
}

impl<Aind: AbsInd> StructureContract for SymbolicTensor<Aind> {
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

// impl<Const> Shadowable<Const> for SymbolicTensor {}

#[allow(dead_code)]
impl<Aind: AbsInd> SymbolicTensor<Aind> {
    pub fn from_named<N>(structure: &N) -> Option<Self>
    where
        N: ToSymbolic + HasName + TensorStructure<Slot = LibrarySlot<Aind>>,
        N::Name: IntoSymbol + Clone,
        N::Args: IntoArgs,
    {
        let permuted_structure = PermutedStructure::from(structure.external_structure());
        let is_metric = structure.name()?.ref_into_symbol() == ETS.metric;
        Some(SymbolicTensor {
            expression: structure.to_symbolic(Some(permuted_structure.index_permutation))?,
            is_metric,
            is_composite: false,
            structure: permuted_structure.structure,
        })
    }

    pub fn from_permuted<N>(structure: &PermutedStructure<N>) -> Option<Self>
    where
        N: ToSymbolic + HasName + TensorStructure<Slot = LibrarySlot<Aind>>,
        N::Name: IntoSymbol + Clone,
        N::Args: IntoArgs,
    {
        let permuted_structure = PermutedStructure::from(structure.structure.external_structure());
        let is_metric = structure.structure.name()?.ref_into_symbol() == ETS.metric;
        Some(SymbolicTensor {
            expression: structure
                .structure
                .to_symbolic(Some(structure.index_permutation.clone()))?,
            is_composite: false,
            is_metric,
            structure: permuted_structure.structure,
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
        Network::<
            NetworkStore<MixedTensor<f64, ShadowedStructure<AbstractIndex>>, Atom>,
            ExplicitKey<AbstractIndex>,
            Symbol,
            AbstractIndex,
        >::try_from_view(
            self.expression.as_view(),
            library,
            &ParseSettings::default(),
        )
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

impl<Aind: AbsInd> HasName for SymbolicTensor<Aind> {
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
impl<Aind: AbsInd> Contract<SymbolicTensor<Aind>> for SymbolicTensor<Aind> {
    type LCM = SymbolicTensor<Aind>;
    fn contract(&self, other: &SymbolicTensor<Aind>) -> Result<Self::LCM, ContractionError> {
        let (out, _, _, _) = self.merge(other)?;
        Ok(out)
    }
}

pub type SymbolicNet<Aind> =
    Network<NetworkStore<SymbolicTensor<Aind>, Atom>, DummyKey, Symbol, Aind>;

// pub type HepTensor<Aind> = MixedTensor<f64, ShadowedStructure<Aind>>;

// pub type ParamNet<Aind> =
//     Network<NetworkStore<ParamTensor<SymbolicTensor<Aind>, Atom>, DummyKey, Symbol, Aind>;
pub trait SymbolicNetExt<Aind: AbsInd + DummyAind + ParseableAind + 'static> {
    fn snapshot_dot(&self) -> String;
    fn simple_execute<CStrat>(self) -> Atom
    where
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, CStrat, LCM = SymbolicTensor<Aind>>;
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

    fn simple_execute<CStrat>(mut self) -> Atom
    where
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, CStrat, LCM = SymbolicTensor<Aind>>,
    {
        let lib = DummyLibrary::<SymbolicTensor<Aind>>::new();

        self.execute::<Sequential, SmallestDegree<CStrat>, _, _, _>(&lib, &Wrap {})
            .unwrap();

        match self.result().unwrap() {
            ExecutionResult::One => Atom::num(1),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(a) => match a {
                TensorOrScalarOrKey::Key { .. } => panic!("aaa"),
                TensorOrScalarOrKey::Scalar(s) => s.clone(),
                TensorOrScalarOrKey::Tensor { tensor, .. } => tensor.expression.clone(),
            },
        }
    }
}

pub trait SymbolicNetParse {
    #[allow(clippy::result_large_err)]
    fn parse_to_symbolic_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>>;
}

impl SymbolicNetParse for Atom {
    fn parse_to_symbolic_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        self.as_view().parse_to_symbolic_net::<Aind>(settings)
    }
}

impl SymbolicNetParse for AtomView<'_> {
    fn parse_to_symbolic_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        let lib = DummyLibrary::<SymbolicTensor<Aind>>::new();

        SymbolicNet::<Aind>::try_from_view::<SymbolicTensor<Aind>, _>(*self, &lib, settings)
    }
}

pub mod contract;

impl<Aind: AbsInd + ParseableAind> Concretize<SymbolicTensor<Aind>> for ShadowedStructure<Aind> {
    fn concretize(self, perm: Option<Permutation>) -> SymbolicTensor<Aind> {
        let is_metric = self.name().unwrap() == ETS.metric;
        SymbolicTensor {
            expression: self.to_symbolic(perm).unwrap(),
            is_composite: false,
            is_metric,
            structure: self.structure,
        }
    }
}

impl<Aind: AbsInd + ParseableAind> Concretize<SymbolicTensor<Aind>>
    for TensorShell<ShadowedStructure<Aind>>
{
    fn concretize(self, perm: Option<Permutation>) -> SymbolicTensor<Aind> {
        let is_metric = self.name().unwrap() == ETS.metric;
        SymbolicTensor {
            expression: self.to_symbolic(perm).unwrap(),
            is_composite: false,
            is_metric,
            structure: self.structure.structure,
        }
    }
}

impl<Aind: AbsInd + ParseableAind> Concretize<SymbolicTensor<Aind>>
    for TensorShell<SymbolicTensor<Aind>>
{
    fn concretize(self, perm: Option<Permutation>) -> SymbolicTensor<Aind> {
        if let Some(p) = perm {
            if p.is_identity() {
                self.structure
            } else {
                panic!("Cannot concretize with a permuation")
            }
        } else {
            self.structure
        }
    }
}

impl<Aind: AbsInd> std::ops::Neg for SymbolicTensor<Aind> {
    type Output = SymbolicTensor<Aind>;

    fn neg(self) -> Self::Output {
        Self {
            expression: -self.expression,
            is_composite: true,
            is_metric: self.is_metric,
            structure: self.structure,
        }
    }
}

impl<Aind: AbsInd> AddAssign<SymbolicTensor<Aind>> for SymbolicTensor<Aind> {
    fn add_assign(&mut self, rhs: SymbolicTensor<Aind>) {
        debug_assert_eq!(self.structure, rhs.structure);
        self.expression += rhs.expression;
        self.is_composite = true;
        self.is_metric = false;
    }
}

impl<Aind: AbsInd> AddAssign<&SymbolicTensor<Aind>> for SymbolicTensor<Aind> {
    fn add_assign(&mut self, rhs: &SymbolicTensor<Aind>) {
        debug_assert_eq!(self.structure, rhs.structure);
        self.expression += &rhs.expression;
        self.is_composite = true;
        self.is_metric = false;
    }
}

impl<Aind: AbsInd> ScalarMul<Atom> for SymbolicTensor<Aind> {
    type Output = SymbolicTensor<Aind>;

    fn scalar_mul(&self, rhs: &Atom) -> Option<Self::Output> {
        Some(Self {
            expression: rhs * &self.expression,
            is_composite: true,
            is_metric: false,
            structure: self.structure.clone(),
        })
    }
}

impl<Aind: AbsInd> std::fmt::Display for SymbolicTensor<Aind> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.expression)
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    use spenso::network::library::symbolic::ETS;
    use symbolica::parse;

    #[test]
    fn parse() {
        let _ = ETS.metric;
        let expr = parse!("g(mink(4,6),mink(4,7))");

        let structure = SymbolicTensor::from_permuted(
            &PermutedStructure::<ShadowedStructure<AbstractIndex>>::try_from(expr).unwrap(),
        )
        .unwrap();

        Contract::<SymbolicTensor, ()>::contract(&structure, &structure).unwrap();
    }
}
