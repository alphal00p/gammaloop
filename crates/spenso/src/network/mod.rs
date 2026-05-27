use graph::{
    NAdd, NMul, NetworkEdge, NetworkGraph, NetworkLeaf, NetworkNode, NetworkOp, NetworkOperation,
    TensorTerm,
};
use linnet::half_edge::{
    NodeIndex,
    subgraph::{Inclusion, SuBitGraph, SubSetLike},
};
use profile::{Counter, Timer};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use library::{Library, LibraryError};

use crate::algebra::{ScalarMul, algebraic_traits::RefOne};
use crate::contraction::{Contract, Trace};
#[cfg(feature = "shadowing")]
use crate::iterators::IteratableTensor;
use crate::network::library::{DummyKey, FunctionLibrary, FunctionLibraryError, LibraryTensor};
use crate::structure::abstract_index::AbstractIndex;
#[cfg(feature = "shadowing")]
use crate::structure::concrete_index::{ConcreteIndex, FlatIndex};
use crate::structure::permuted::PermuteTensor;
// use crate::shadowing::Concretize;
#[cfg(feature = "shadowing")]
use crate::structure::StructureContract;
use crate::structure::representation::LibrarySlot;
use crate::structure::slot::{AbsInd, IsAbstractSlot};
use crate::structure::{HasName, PermutedStructure, StructureError, TensorShell};
use std::borrow::Cow;
#[cfg(feature = "shadowing")]
use std::collections::HashMap;
use std::fmt::Display;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};
use store::{
    NetworkStore, NetworkStoreAccess, NetworkStoreOverlay, TensorScalarStore,
    TensorScalarStoreMapping,
};
use thiserror::Error;
// use log::trace;

#[cfg(feature = "shadowing")]
use crate::algebra::complex::RealOrComplexRef;
#[cfg(feature = "shadowing")]
use crate::tensors::parametric::{
    AtomViewOrConcrete, ConcreteOrParam, ParamOrConcrete, ParamTensor,
};
use crate::tensors::{
    complex::RealOrComplexTensor,
    data::{DataTensor, DenseTensor, SparseTensor},
};
use crate::{
    contraction::ContractionError,
    structure::{CastStructure, HasStructure, ScalarTensor, TensorStructure},
};
use eyre::eyre;

#[cfg(feature = "shadowing")]
use symbolica::atom::{Atom, AtomCore, AtomView};
#[cfg(feature = "shadowing")]
use symbolica::id::AliasedAtom;

#[cfg(feature = "shadowing")]
pub mod tags;
#[cfg(feature = "shadowing")]
use tags::scalar_store_alias;
// use eyre::Result;

use std::{convert::Infallible, fmt::Debug};

const LARGE_SUM_PROFILE_THRESHOLD: usize = 32;
const MIN_LAZY_TENSOR_SUM_TERMS: usize = 2;
#[cfg(feature = "shadowing")]
const MIN_LAZY_FUSED_NUMERIC_CONTRACT_TERMS: usize = 8;

pub struct FastTensorSumContractTerm<T, Sc> {
    pub tensor: T,
    pub scalar: Option<Sc>,
}

pub enum FastTensorSumContract<T, Sc> {
    Materialized(T),
    Terms(Vec<T>),
    ScaledTerms(Vec<FastTensorSumContractTerm<T, Sc>>),
}

impl<T, Sc> FastTensorSumContract<T, Sc> {
    #[cfg(feature = "shadowing")]
    fn map<U>(self, mut f: impl FnMut(T) -> U) -> FastTensorSumContract<U, Sc> {
        match self {
            FastTensorSumContract::Materialized(tensor) => {
                FastTensorSumContract::Materialized(f(tensor))
            }
            FastTensorSumContract::Terms(terms) => {
                FastTensorSumContract::Terms(terms.into_iter().map(f).collect())
            }
            FastTensorSumContract::ScaledTerms(terms) => FastTensorSumContract::ScaledTerms(
                terms
                    .into_iter()
                    .map(|term| FastTensorSumContractTerm {
                        tensor: f(term.tensor),
                        scalar: term.scalar,
                    })
                    .collect(),
            ),
        }
    }
}

pub trait FastTensorSum: Sized {
    fn fast_tensor_sum(_terms: &[&Self], _sum_start: Option<std::time::Instant>) -> Option<Self> {
        None
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TensorContractionProfile {
    pub entries: usize,
    pub total_terms: usize,
    pub max_terms: usize,
    pub total_bytes: usize,
    pub max_bytes: usize,
    pub common_factor_count: usize,
    pub simple_tensor: bool,
}

impl Default for TensorContractionProfile {
    fn default() -> Self {
        Self::unit()
    }
}

impl TensorContractionProfile {
    pub fn unit() -> Self {
        Self {
            entries: 1,
            total_terms: 1,
            max_terms: 1,
            total_bytes: 1,
            max_bytes: 1,
            common_factor_count: 0,
            simple_tensor: false,
        }
    }

    #[cfg(feature = "shadowing")]
    fn from_param_tensor<S>(tensor: &ParamTensor<S>) -> Self
    where
        S: TensorStructure + Clone,
    {
        let mut profile = Self {
            entries: tensor.tensor.actual_size().max(1),
            total_terms: 0,
            max_terms: 0,
            total_bytes: 0,
            max_bytes: 0,
            common_factor_count: 0,
            simple_tensor: false,
        };
        let mut nonzero_entries = 0usize;
        let mut numeric_entries = 0usize;

        let mut observe = |profile: &mut Self, atom: AtomView<'_>| {
            if !atom.is_zero() {
                nonzero_entries += 1;
                if matches!(atom, AtomView::Num(_)) {
                    numeric_entries += 1;
                }
            }
            let terms = atom.nterms();
            let bytes = atom.get_byte_size();
            profile.total_terms += terms;
            profile.max_terms = profile.max_terms.max(terms);
            profile.total_bytes += bytes;
            profile.max_bytes = profile.max_bytes.max(bytes);
        };

        match &tensor.tensor {
            DataTensor::Dense(dense) => {
                for atom in &dense.data {
                    observe(&mut profile, atom.as_view());
                }
            }
            DataTensor::Sparse(sparse) => {
                for atom in sparse.elements.values() {
                    observe(&mut profile, atom.as_view());
                }
                profile.common_factor_count = common_sparse_atom_factors(sparse).len();
            }
        }

        profile.total_terms = profile.total_terms.max(1);
        profile.max_terms = profile.max_terms.max(1);
        profile.total_bytes = profile.total_bytes.max(1);
        profile.max_bytes = profile.max_bytes.max(1);
        profile.simple_tensor = nonzero_entries > 0 && nonzero_entries == numeric_entries;
        profile
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TensorContractionPairEstimate {
    pub estimated_output_entries: u128,
    pub output_dense_size: u128,
    pub max_output_entry_products: u128,
    pub simple_tensor_penalty: u128,
    pub common_factor_penalty: u128,
}

impl TensorContractionPairEstimate {
    pub fn from_profiles(
        left: TensorContractionProfile,
        right: TensorContractionProfile,
        output_dense_size: u128,
    ) -> Self {
        let entry_work = (left.entries.max(1) as u128).saturating_mul(right.entries.max(1) as u128);
        let output_dense_size = output_dense_size.max(1);
        let estimated_output_entries = output_dense_size.min(entry_work).max(1);
        Self {
            estimated_output_entries,
            output_dense_size,
            max_output_entry_products: div_ceil_u128(entry_work, estimated_output_entries).max(1),
            simple_tensor_penalty: u128::from(!(left.simple_tensor || right.simple_tensor)),
            common_factor_penalty: u128::from(
                left.common_factor_count == 0 && right.common_factor_count == 0,
            ),
        }
    }
}

fn div_ceil_u128(numerator: u128, denominator: u128) -> u128 {
    if denominator == 0 {
        return numerator;
    }
    numerator.div_ceil(denominator)
}

pub trait FastTensorSumContractible<Sc>: Sized {
    fn fast_tensor_sum_contract<CStrat>(
        _terms: &[&Self],
        _other: &Self,
        _terms_on_left: bool,
    ) -> Option<Result<FastTensorSumContract<Self, Sc>, ContractionError>> {
        None
    }

    fn contraction_profile(&self) -> TensorContractionProfile {
        TensorContractionProfile::unit()
    }

    fn contraction_pair_estimate(
        &self,
        other: &Self,
        _left_matches: &[bool],
        _right_matches: &[bool],
        output_dense_size: u128,
        _exact_join_limit: usize,
    ) -> TensorContractionPairEstimate {
        TensorContractionPairEstimate::from_profiles(
            self.contraction_profile(),
            other.contraction_profile(),
            output_dense_size,
        )
    }
}

pub trait TensorCommonFactor<Sc>: Sized {
    fn split_common_factor(&self) -> Option<(Self, Sc)> {
        None
    }
}

impl<T, S> FastTensorSum for DenseTensor<T, S> {}
impl<T, S, Sc> FastTensorSumContractible<Sc> for DenseTensor<T, S> {}
impl<T, S, Sc> TensorCommonFactor<Sc> for DenseTensor<T, S> {}
impl<T, S> FastTensorSum for SparseTensor<T, S> {}
impl<T, S, Sc> FastTensorSumContractible<Sc> for SparseTensor<T, S> {}
impl<T, S, Sc> TensorCommonFactor<Sc> for SparseTensor<T, S> {}
impl<T, S> FastTensorSum for DataTensor<T, S> {}
impl<T, S, Sc> FastTensorSumContractible<Sc> for DataTensor<T, S> {}
impl<T, S, Sc> TensorCommonFactor<Sc> for DataTensor<T, S> {}
impl<T, S: TensorStructure> FastTensorSum for RealOrComplexTensor<T, S> {}
impl<T, S: TensorStructure, Sc> FastTensorSumContractible<Sc> for RealOrComplexTensor<T, S> {}
impl<T, S: TensorStructure, Sc> TensorCommonFactor<Sc> for RealOrComplexTensor<T, S> {}
impl<S> FastTensorSum for TensorShell<S> {}
impl<S, Sc> FastTensorSumContractible<Sc> for TensorShell<S> {}
impl<S, Sc> TensorCommonFactor<Sc> for TensorShell<S> {}

#[cfg(feature = "shadowing")]
impl<S> FastTensorSum for ParamTensor<S>
where
    S: TensorStructure + Clone + Sync,
{
    fn fast_tensor_sum(terms: &[&Self], sum_start: Option<std::time::Instant>) -> Option<Self> {
        if terms.len() < 2 {
            return None;
        }

        let first = *terms.first()?;
        let structure = first.tensor.structure().clone();
        let mut entries = HashMap::<FlatIndex, Vec<AtomView<'_>>>::new();
        let mut input_entries = 0usize;

        for term in terms {
            let DataTensor::Sparse(sparse) = &term.tensor else {
                return None;
            };
            if !sparse.zero.as_view().is_zero() {
                return None;
            }

            input_entries += sparse.elements.len();
            for (index, atom) in &sparse.elements {
                if !atom.as_view().is_zero() {
                    entries.entry(*index).or_default().push(atom.as_view());
                }
            }
        }

        let log_fast_sum = profile::verbose() || sum_start.is_some();
        let fast_start = log_fast_sum.then(std::time::Instant::now);
        if log_fast_sum && let Some(sum_start) = sum_start {
            eprintln!(
                "spenso_profile execute.fast_tensor_sum_start terms={} input_entries={} grouped_entries={} elapsed_ms={:.3}",
                terms.len(),
                input_entries,
                entries.len(),
                sum_start.elapsed().as_secs_f64() * 1000.0,
            );
        }

        let merged = entries
            .into_iter()
            .collect::<Vec<_>>()
            .into_par_iter()
            .filter_map(|(index, atoms)| {
                let atom = match atoms.len() {
                    0 => Atom::Zero,
                    1 => atoms
                        .into_iter()
                        .next()
                        .expect("single atom exists")
                        .to_owned(),
                    _ => Atom::add_many(&atoms),
                };
                (!atom.as_view().is_zero()).then_some((index, atom))
            })
            .collect::<Vec<_>>();

        let mut elements = HashMap::with_capacity(merged.len());
        for (index, atom) in merged {
            elements.insert(index, atom);
        }

        if let Some(start) = fast_start {
            eprintln!(
                "spenso_profile execute.fast_tensor_sum_done terms={} input_entries={} output_entries={} elapsed_ms={:.3}",
                terms.len(),
                input_entries,
                elements.len(),
                start.elapsed().as_secs_f64() * 1000.0,
            );
        }

        Some(ParamTensor::composite(DataTensor::Sparse(SparseTensor {
            elements,
            zero: Atom::Zero,
            structure,
        })))
    }
}

#[cfg(feature = "shadowing")]
impl<S> FastTensorSumContractible<Atom> for ParamTensor<S>
where
    S: TensorStructure + Clone + Send + Sync + StructureContract,
    S::Slot: Send + Sync,
{
    fn contraction_profile(&self) -> TensorContractionProfile {
        TensorContractionProfile::from_param_tensor(self)
    }

    fn contraction_pair_estimate(
        &self,
        other: &Self,
        left_matches: &[bool],
        right_matches: &[bool],
        output_dense_size: u128,
        exact_join_limit: usize,
    ) -> TensorContractionPairEstimate {
        sparse_atom_pair_estimate(
            self,
            other,
            left_matches,
            right_matches,
            output_dense_size,
            exact_join_limit,
        )
    }

    fn fast_tensor_sum_contract<CStrat>(
        terms: &[&Self],
        other: &Self,
        terms_on_left: bool,
    ) -> Option<Result<FastTensorSumContract<Self, Atom>, ContractionError>> {
        fused_numeric_tensor_sum_contract(terms, other, terms_on_left)
    }
}

#[cfg(feature = "shadowing")]
impl<S> TensorCommonFactor<Atom> for ParamTensor<S>
where
    S: TensorStructure + Clone,
{
    fn split_common_factor(&self) -> Option<(Self, Atom)> {
        let DataTensor::Sparse(sparse) = &self.tensor else {
            return None;
        };
        if !sparse.zero.as_view().is_zero() {
            return None;
        }

        let common_factors = common_sparse_atom_factors(sparse);
        if common_factors.is_empty() {
            return None;
        }

        let factor = atom_product(common_factors.iter().cloned());
        if factor == Atom::num(1) {
            return None;
        }

        let elements = sparse
            .elements
            .iter()
            .filter_map(|(index, atom)| {
                let reduced = atom_without_factors(atom.as_view(), &common_factors);
                (!reduced.as_view().is_zero()).then_some((*index, reduced))
            })
            .collect::<HashMap<_, _>>();

        Some((
            ParamTensor::composite(DataTensor::Sparse(SparseTensor {
                elements,
                zero: Atom::Zero,
                structure: sparse.structure.clone(),
            })),
            factor,
        ))
    }
}

#[cfg(feature = "shadowing")]
fn sparse_atom_pair_estimate<S>(
    left: &ParamTensor<S>,
    right: &ParamTensor<S>,
    left_matches: &[bool],
    right_matches: &[bool],
    output_dense_size: u128,
    exact_join_limit: usize,
) -> TensorContractionPairEstimate
where
    S: TensorStructure + Clone,
{
    let left_profile = TensorContractionProfile::from_param_tensor(left);
    let right_profile = TensorContractionProfile::from_param_tensor(right);
    let fallback = TensorContractionPairEstimate::from_profiles(
        left_profile,
        right_profile,
        output_dense_size,
    );

    let DataTensor::Sparse(left_sparse) = &left.tensor else {
        return fallback;
    };
    let DataTensor::Sparse(right_sparse) = &right.tensor else {
        return fallback;
    };
    if !left_sparse.zero.as_view().is_zero() || !right_sparse.zero.as_view().is_zero() {
        return fallback;
    }

    let left_groups = sparse_free_keys_by_match(left_sparse, left_matches);
    let right_groups = sparse_free_keys_by_match(right_sparse, right_matches);
    if left_groups.is_empty() || right_groups.is_empty() {
        return TensorContractionPairEstimate {
            estimated_output_entries: 1,
            output_dense_size: output_dense_size.max(1),
            max_output_entry_products: 1,
            simple_tensor_penalty: u128::from(
                !(left_profile.simple_tensor || right_profile.simple_tensor),
            ),
            common_factor_penalty: u128::from(
                left_profile.common_factor_count == 0 && right_profile.common_factor_count == 0,
            ),
        };
    }

    let mut join_cardinality = 0u128;
    let mut max_key_products = 0u128;
    let mut exact_combinations = 0usize;
    for (key, left_free_keys) in &left_groups {
        let Some(right_free_keys) = right_groups.get(key) else {
            continue;
        };
        let key_products =
            (left_free_keys.len() as u128).saturating_mul(right_free_keys.len() as u128);
        join_cardinality = join_cardinality.saturating_add(key_products);
        max_key_products = max_key_products.max(key_products);
        exact_combinations = exact_combinations
            .saturating_add(left_free_keys.len().saturating_mul(right_free_keys.len()));
    }

    if join_cardinality == 0 {
        return TensorContractionPairEstimate {
            estimated_output_entries: 1,
            output_dense_size: output_dense_size.max(1),
            max_output_entry_products: 1,
            simple_tensor_penalty: u128::from(
                !(left_profile.simple_tensor || right_profile.simple_tensor),
            ),
            common_factor_penalty: u128::from(
                left_profile.common_factor_count == 0 && right_profile.common_factor_count == 0,
            ),
        };
    }

    let output_dense_size = output_dense_size.max(1);
    let (estimated_output_entries, max_output_entry_products) = if exact_combinations
        <= exact_join_limit
    {
        let mut output_counts = HashMap::<(Vec<ConcreteIndex>, Vec<ConcreteIndex>), u128>::new();
        for (key, left_free_keys) in &left_groups {
            let Some(right_free_keys) = right_groups.get(key) else {
                continue;
            };
            for left_key in left_free_keys {
                for right_key in right_free_keys {
                    *output_counts
                        .entry((left_key.clone(), right_key.clone()))
                        .or_default() += 1;
                }
            }
        }
        let estimated_output_entries = (output_counts.len() as u128).clamp(1, output_dense_size);
        let max_output_entry_products = output_counts.values().copied().max().unwrap_or(1);
        (estimated_output_entries, max_output_entry_products)
    } else {
        let estimated_output_entries = output_dense_size.min(join_cardinality).max(1);
        let max_output_entry_products = max_key_products
            .max(div_ceil_u128(join_cardinality, estimated_output_entries))
            .max(1);
        (estimated_output_entries, max_output_entry_products)
    };

    TensorContractionPairEstimate {
        estimated_output_entries,
        output_dense_size,
        max_output_entry_products,
        simple_tensor_penalty: u128::from(
            !(left_profile.simple_tensor || right_profile.simple_tensor),
        ),
        common_factor_penalty: u128::from(
            left_profile.common_factor_count == 0 && right_profile.common_factor_count == 0,
        ),
    }
}

#[cfg(feature = "shadowing")]
fn sparse_free_keys_by_match<S>(
    sparse: &SparseTensor<Atom, S>,
    matches: &[bool],
) -> HashMap<Vec<ConcreteIndex>, Vec<Vec<ConcreteIndex>>>
where
    S: TensorStructure,
{
    let mut groups = HashMap::<Vec<ConcreteIndex>, Vec<Vec<ConcreteIndex>>>::new();
    for (expanded, atom) in sparse.iter_expanded() {
        if atom.as_view().is_zero() {
            continue;
        }
        let mut matched_key = Vec::new();
        let mut free_key = Vec::new();
        for (axis, index) in expanded.indices.iter().copied().enumerate() {
            if matches.get(axis).copied().unwrap_or(false) {
                matched_key.push(index);
            } else {
                free_key.push(index);
            }
        }
        groups.entry(matched_key).or_default().push(free_key);
    }
    groups
}

#[cfg(feature = "shadowing")]
impl<C, S> FastTensorSum for ParamOrConcrete<C, S>
where
    S: TensorStructure + Clone + Sync,
{
    fn fast_tensor_sum(terms: &[&Self], sum_start: Option<std::time::Instant>) -> Option<Self> {
        let param_terms = terms
            .iter()
            .map(|term| match term {
                ParamOrConcrete::Param(param) => Some(param),
                ParamOrConcrete::Concrete(_) => None,
            })
            .collect::<Option<Vec<_>>>()?;

        <ParamTensor<S> as FastTensorSum>::fast_tensor_sum(&param_terms, sum_start)
            .map(ParamOrConcrete::Param)
    }
}

#[cfg(feature = "shadowing")]
impl<C, S> FastTensorSumContractible<Atom> for ParamOrConcrete<C, S>
where
    S: TensorStructure + Clone + Send + Sync + StructureContract,
    S::Slot: Send + Sync,
{
    fn contraction_profile(&self) -> TensorContractionProfile {
        match self {
            ParamOrConcrete::Param(param) => TensorContractionProfile::from_param_tensor(param),
            ParamOrConcrete::Concrete(_) => TensorContractionProfile::unit(),
        }
    }

    fn contraction_pair_estimate(
        &self,
        other: &Self,
        left_matches: &[bool],
        right_matches: &[bool],
        output_dense_size: u128,
        exact_join_limit: usize,
    ) -> TensorContractionPairEstimate {
        match (self, other) {
            (ParamOrConcrete::Param(left), ParamOrConcrete::Param(right)) => {
                sparse_atom_pair_estimate(
                    left,
                    right,
                    left_matches,
                    right_matches,
                    output_dense_size,
                    exact_join_limit,
                )
            }
            _ => TensorContractionPairEstimate::from_profiles(
                <Self as FastTensorSumContractible<Atom>>::contraction_profile(self),
                <Self as FastTensorSumContractible<Atom>>::contraction_profile(other),
                output_dense_size,
            ),
        }
    }

    fn fast_tensor_sum_contract<CStrat>(
        terms: &[&Self],
        other: &Self,
        terms_on_left: bool,
    ) -> Option<Result<FastTensorSumContract<Self, Atom>, ContractionError>> {
        let param_terms = terms
            .iter()
            .map(|term| match term {
                ParamOrConcrete::Param(param) => Some(param),
                ParamOrConcrete::Concrete(_) => None,
            })
            .collect::<Option<Vec<_>>>()?;
        let ParamOrConcrete::Param(other) = other else {
            return None;
        };

        <ParamTensor<S> as FastTensorSumContractible<Atom>>::fast_tensor_sum_contract::<CStrat>(
            &param_terms,
            other,
            terms_on_left,
        )
        .map(|result| result.map(|contract| contract.map(ParamOrConcrete::Param)))
    }
}

#[cfg(feature = "shadowing")]
impl<C, S, CSc> FastTensorSumContractible<ConcreteOrParam<CSc>> for ParamOrConcrete<C, S>
where
    S: TensorStructure + Clone + Send + Sync + StructureContract,
    S::Slot: Send + Sync,
{
    fn contraction_profile(&self) -> TensorContractionProfile {
        match self {
            ParamOrConcrete::Param(param) => TensorContractionProfile::from_param_tensor(param),
            ParamOrConcrete::Concrete(_) => TensorContractionProfile::unit(),
        }
    }

    fn fast_tensor_sum_contract<CStrat>(
        terms: &[&Self],
        other: &Self,
        terms_on_left: bool,
    ) -> Option<Result<FastTensorSumContract<Self, ConcreteOrParam<CSc>>, ContractionError>> {
        let param_terms = terms
            .iter()
            .map(|term| match term {
                ParamOrConcrete::Param(param) => Some(param),
                ParamOrConcrete::Concrete(_) => None,
            })
            .collect::<Option<Vec<_>>>()?;
        let ParamOrConcrete::Param(other) = other else {
            return None;
        };

        <ParamTensor<S> as FastTensorSumContractible<Atom>>::fast_tensor_sum_contract::<CStrat>(
            &param_terms,
            other,
            terms_on_left,
        )
        .map(|result| {
            result.map(|contract| match contract {
                FastTensorSumContract::Materialized(tensor) => {
                    FastTensorSumContract::Materialized(ParamOrConcrete::Param(tensor))
                }
                FastTensorSumContract::Terms(terms) => FastTensorSumContract::Terms(
                    terms.into_iter().map(ParamOrConcrete::Param).collect(),
                ),
                FastTensorSumContract::ScaledTerms(terms) => FastTensorSumContract::ScaledTerms(
                    terms
                        .into_iter()
                        .map(|term| FastTensorSumContractTerm {
                            tensor: ParamOrConcrete::Param(term.tensor),
                            scalar: term.scalar.map(ConcreteOrParam::Param),
                        })
                        .collect(),
                ),
            })
        })
    }
}

#[cfg(feature = "shadowing")]
impl<C, S> TensorCommonFactor<Atom> for ParamOrConcrete<C, S>
where
    S: TensorStructure + Clone,
{
    fn split_common_factor(&self) -> Option<(Self, Atom)> {
        match self {
            ParamOrConcrete::Param(param) => param
                .split_common_factor()
                .map(|(tensor, factor)| (ParamOrConcrete::Param(tensor), factor)),
            ParamOrConcrete::Concrete(_) => None,
        }
    }
}

#[cfg(feature = "shadowing")]
impl<C, S, CSc> TensorCommonFactor<ConcreteOrParam<CSc>> for ParamOrConcrete<C, S>
where
    S: TensorStructure + Clone,
{
    fn split_common_factor(&self) -> Option<(Self, ConcreteOrParam<CSc>)> {
        match self {
            ParamOrConcrete::Param(param) => param.split_common_factor().map(|(tensor, factor)| {
                (
                    ParamOrConcrete::Param(tensor),
                    ConcreteOrParam::Param(factor),
                )
            }),
            ParamOrConcrete::Concrete(_) => None,
        }
    }
}

#[cfg(feature = "shadowing")]
struct NumericContractEntry {
    index: Vec<ConcreteIndex>,
    coefficient: Atom,
    coefficient_class: usize,
    left_metric_negative: bool,
}

#[cfg(feature = "shadowing")]
struct FusedNumericContractPlan<S>
where
    S: TensorStructure,
{
    terms_on_left: bool,
    result_structure: S,
    sum_structure: S,
    numeric_entries: Vec<NumericContractEntry>,
    coefficient_classes: Vec<Atom>,
    numeric_by_key: HashMap<Vec<ConcreteIndex>, Vec<usize>>,
    sum_match_axes: Vec<usize>,
    left_slots: Vec<S::Slot>,
    right_slots: Vec<S::Slot>,
    left_match_flags: Vec<bool>,
    right_match_flags: Vec<bool>,
}

#[cfg(feature = "shadowing")]
fn fused_numeric_tensor_sum_contract<S>(
    terms: &[&ParamTensor<S>],
    numeric: &ParamTensor<S>,
    terms_on_left: bool,
) -> Option<Result<FastTensorSumContract<ParamTensor<S>, Atom>, ContractionError>>
where
    S: TensorStructure + Clone + Send + Sync + StructureContract,
    S::Slot: Send + Sync,
{
    if terms.is_empty() {
        return None;
    }

    let first = *terms.first()?;
    let DataTensor::Sparse(first_sparse) = &first.tensor else {
        return None;
    };
    if !first_sparse.zero.as_view().is_zero() {
        return None;
    }
    let DataTensor::Sparse(numeric_sparse) = &numeric.tensor else {
        return None;
    };
    if !numeric_sparse.zero.as_view().is_zero()
        || !numeric_sparse
            .elements
            .values()
            .all(|atom| matches!(atom.as_view(), AtomView::Num(_)))
    {
        return None;
    }

    let sum_structure = first_sparse.structure.clone();
    let sum_slots = sum_structure.external_structure();
    if terms
        .iter()
        .any(|term| term.tensor.structure().external_structure() != sum_slots)
    {
        return None;
    }

    let numeric_structure = numeric_sparse.structure.clone();
    let (left_structure, right_structure) = if terms_on_left {
        (&sum_structure, &numeric_structure)
    } else {
        (&numeric_structure, &sum_structure)
    };
    if let Ok((_, left_matches, _, _)) = left_structure.merge(right_structure)
        && left_matches.n_included() == 0
    {
        return None;
    }

    Some((|| {
        let plan = FusedNumericContractPlan::new(
            first_sparse.structure.clone(),
            numeric_sparse,
            terms_on_left,
        )?;

        if terms.len() == 1 {
            plan.contract_term(terms[0])
                .map(|term| FastTensorSumContract::ScaledTerms(vec![term]))
        } else if terms.len() >= MIN_LAZY_FUSED_NUMERIC_CONTRACT_TERMS {
            let start = profile::enabled().then(std::time::Instant::now);
            let contracted_terms = terms
                .par_iter()
                .map(|term| plan.contract_term(term))
                .collect::<Result<Vec<_>, _>>()?;
            if let Some(start) = start {
                let output_entries = contracted_terms
                    .iter()
                    .map(|term| term.tensor.tensor.actual_size())
                    .sum::<usize>();
                let scaled_terms = contracted_terms
                    .iter()
                    .filter(|term| term.scalar.is_some())
                    .count();
                eprintln!(
                    "spenso_profile product.fused_numeric_tensor_sum_lazy_contract terms={} numeric_entries={} output_terms={} scaled_terms={} output_entries={} elapsed_ms={:.3}",
                    terms.len(),
                    plan.numeric_entries.len(),
                    contracted_terms.len(),
                    scaled_terms,
                    output_entries,
                    start.elapsed().as_secs_f64() * 1000.0,
                );
            }
            Ok(FastTensorSumContract::ScaledTerms(contracted_terms))
        } else {
            fused_numeric_tensor_sum_contract_impl(terms, plan)
                .map(FastTensorSumContract::Materialized)
        }
    })())
}

#[cfg(feature = "shadowing")]
fn fused_numeric_tensor_sum_contract_impl<S>(
    terms: &[&ParamTensor<S>],
    plan: FusedNumericContractPlan<S>,
) -> Result<ParamTensor<S>, ContractionError>
where
    S: TensorStructure + Clone + Send + Sync + StructureContract,
    S::Slot: Send + Sync,
{
    let mut groups = HashMap::<(FlatIndex, usize, bool), Vec<AtomView<'_>>>::new();

    for term in terms {
        let DataTensor::Sparse(sum_sparse) = &term.tensor else {
            return Err(eyre!("fused tensor sum term is not sparse").into());
        };
        plan.collect_borrowed_groups(sum_sparse, &mut groups)?;
    }

    if profile::enabled() {
        eprintln!(
            "spenso_profile product.fused_numeric_tensor_sum_contract terms={} numeric_entries={} coefficient_classes={} groups={}",
            terms.len(),
            plan.numeric_entries.len(),
            plan.coefficient_classes.len(),
            groups.len(),
        );
    }

    plan.finish_groups(groups)
}

#[cfg(feature = "shadowing")]
impl<S> FusedNumericContractPlan<S>
where
    S: TensorStructure + Clone + Send + Sync + StructureContract,
    S::Slot: Send + Sync,
{
    fn new(
        sum_structure: S,
        numeric: &SparseTensor<Atom, S>,
        terms_on_left: bool,
    ) -> Result<Self, ContractionError> {
        let numeric_structure = numeric.structure.clone();
        let (left_structure, right_structure) = if terms_on_left {
            (&sum_structure, &numeric_structure)
        } else {
            (&numeric_structure, &sum_structure)
        };
        let (result_structure, left_matches, right_matches, _) =
            left_structure.merge(right_structure)?;
        if left_matches.n_included() == 0 {
            return Err(
                eyre!("fused tensor sum contract requires at least one shared slot").into(),
            );
        }

        let left_match_axes = left_matches
            .included_iter()
            .map(|axis| axis.0)
            .collect::<Vec<_>>();
        let right_match_axes = right_matches
            .included_iter()
            .map(|axis| axis.0)
            .collect::<Vec<_>>();
        let (sum_match_axes, numeric_match_axes) = if terms_on_left {
            (left_match_axes.clone(), right_match_axes.clone())
        } else {
            (right_match_axes.clone(), left_match_axes.clone())
        };

        let mut numeric_entries = collect_numeric_contract_entries(
            numeric,
            &numeric_match_axes,
            if terms_on_left {
                None
            } else {
                Some(&left_match_axes)
            },
        )?;
        let mut numeric_by_key = HashMap::<Vec<ConcreteIndex>, Vec<usize>>::new();
        for (entry_index, entry) in numeric_entries.iter().enumerate() {
            numeric_by_key
                .entry(component_key(&entry.index, &numeric_match_axes))
                .or_default()
                .push(entry_index);
        }
        let mut coefficient_classes = Vec::<Atom>::new();
        for entry in &mut numeric_entries {
            let coefficient_class = coefficient_classes
                .iter()
                .position(|coefficient| coefficient == &entry.coefficient)
                .unwrap_or_else(|| {
                    coefficient_classes.push(entry.coefficient.clone());
                    coefficient_classes.len() - 1
                });
            entry.coefficient_class = coefficient_class;
        }

        let left_slots = left_structure.external_structure();
        let right_slots = right_structure.external_structure();
        let left_match_flags = match_flags(left_structure.order(), &left_match_axes);
        let right_match_flags = match_flags(right_structure.order(), &right_match_axes);

        Ok(Self {
            terms_on_left,
            result_structure,
            sum_structure,
            numeric_entries,
            coefficient_classes,
            numeric_by_key,
            sum_match_axes,
            left_slots,
            right_slots,
            left_match_flags,
            right_match_flags,
        })
    }

    fn contract_term(
        &self,
        term: &ParamTensor<S>,
    ) -> Result<FastTensorSumContractTerm<ParamTensor<S>, Atom>, ContractionError> {
        let DataTensor::Sparse(sum_sparse) = &term.tensor else {
            return Err(eyre!("fused tensor sum term is not sparse").into());
        };

        let common_factors = common_sparse_atom_factors(sum_sparse);
        if !common_factors.is_empty() {
            return self.contract_factored_term(sum_sparse, common_factors);
        }

        let mut groups = HashMap::<(FlatIndex, usize, bool), Vec<AtomView<'_>>>::new();
        self.collect_borrowed_groups(sum_sparse, &mut groups)?;

        Ok(FastTensorSumContractTerm {
            tensor: self.finish_groups(groups)?,
            scalar: None,
        })
    }

    fn collect_borrowed_groups<'a>(
        &self,
        sum_sparse: &'a SparseTensor<Atom, S>,
        groups: &mut HashMap<(FlatIndex, usize, bool), Vec<AtomView<'a>>>,
    ) -> Result<(), ContractionError> {
        for (sum_index, atom) in sum_sparse.iter_expanded() {
            if atom.as_view().is_zero() {
                continue;
            }
            let key = component_key(&sum_index.indices, &self.sum_match_axes);
            let Some(matching_numeric_entries) = self.numeric_by_key.get(&key) else {
                continue;
            };
            let sum_left_metric_negative = if self.terms_on_left {
                metric_negative(
                    &self.sum_structure,
                    &self.sum_match_axes,
                    &sum_index.indices,
                )?
            } else {
                false
            };

            for numeric_entry_index in matching_numeric_entries {
                let numeric_entry = &self.numeric_entries[*numeric_entry_index];
                let (left_index, right_index, left_negative) = if self.terms_on_left {
                    (
                        sum_index.indices.as_slice(),
                        numeric_entry.index.as_slice(),
                        sum_left_metric_negative,
                    )
                } else {
                    (
                        numeric_entry.index.as_slice(),
                        sum_index.indices.as_slice(),
                        numeric_entry.left_metric_negative,
                    )
                };
                let flat_index = self.contract_flat_index(left_index, right_index)?;
                groups
                    .entry((flat_index, numeric_entry.coefficient_class, left_negative))
                    .or_default()
                    .push(atom.as_view());
            }
        }

        Ok(())
    }

    fn contract_factored_term(
        &self,
        sum_sparse: &SparseTensor<Atom, S>,
        common_factors: Vec<Atom>,
    ) -> Result<FastTensorSumContractTerm<ParamTensor<S>, Atom>, ContractionError> {
        let factor = atom_product(common_factors.iter().cloned());
        let mut groups = HashMap::<(FlatIndex, usize, bool), Vec<Atom>>::new();
        for (sum_index, atom) in sum_sparse.iter_expanded() {
            if atom.as_view().is_zero() {
                continue;
            }
            let key = component_key(&sum_index.indices, &self.sum_match_axes);
            let Some(matching_numeric_entries) = self.numeric_by_key.get(&key) else {
                continue;
            };
            let reduced_atom = atom_without_factors(atom.as_view(), &common_factors);
            let sum_left_metric_negative = if self.terms_on_left {
                metric_negative(
                    &self.sum_structure,
                    &self.sum_match_axes,
                    &sum_index.indices,
                )?
            } else {
                false
            };

            for numeric_entry_index in matching_numeric_entries {
                let numeric_entry = &self.numeric_entries[*numeric_entry_index];
                let (left_index, right_index, left_negative) = if self.terms_on_left {
                    (
                        sum_index.indices.as_slice(),
                        numeric_entry.index.as_slice(),
                        sum_left_metric_negative,
                    )
                } else {
                    (
                        numeric_entry.index.as_slice(),
                        sum_index.indices.as_slice(),
                        numeric_entry.left_metric_negative,
                    )
                };
                let flat_index = self.contract_flat_index(left_index, right_index)?;
                groups
                    .entry((flat_index, numeric_entry.coefficient_class, left_negative))
                    .or_default()
                    .push(reduced_atom.clone());
            }
        }

        Ok(FastTensorSumContractTerm {
            tensor: self.finish_owned_groups(groups)?,
            scalar: (factor != Atom::num(1)).then_some(factor),
        })
    }

    fn contract_flat_index(
        &self,
        left_index: &[ConcreteIndex],
        right_index: &[ConcreteIndex],
    ) -> Result<FlatIndex, ContractionError> {
        let result_index = contract_result_index(
            &self.result_structure,
            &self.left_slots,
            &self.right_slots,
            &self.left_match_flags,
            &self.right_match_flags,
            left_index,
            right_index,
        )?;
        Ok(self.result_structure.flat_index(&result_index)?)
    }

    fn finish_groups(
        &self,
        groups: HashMap<(FlatIndex, usize, bool), Vec<AtomView<'_>>>,
    ) -> Result<ParamTensor<S>, ContractionError> {
        let contributions = groups
            .into_iter()
            .collect::<Vec<_>>()
            .into_par_iter()
            .filter_map(|((flat_index, coefficient_class, left_negative), atoms)| {
                let mut contribution = match atoms.len() {
                    0 => Atom::Zero,
                    1 => atoms
                        .into_iter()
                        .next()
                        .expect("single grouped atom exists")
                        .to_owned(),
                    _ => Atom::add_many(&atoms),
                };
                if left_negative {
                    contribution = -contribution;
                }
                contribution = multiply_atom_by_numeric_coefficient(
                    contribution,
                    &self.coefficient_classes[coefficient_class],
                );
                (!contribution.as_view().is_zero()).then_some((flat_index, contribution))
            })
            .collect::<Vec<_>>();

        Ok(ParamTensor::composite(DataTensor::Sparse(SparseTensor {
            elements: self.collect_output_groups(contributions),
            zero: Atom::Zero,
            structure: self.result_structure.clone(),
        })))
    }

    fn finish_owned_groups(
        &self,
        groups: HashMap<(FlatIndex, usize, bool), Vec<Atom>>,
    ) -> Result<ParamTensor<S>, ContractionError> {
        let contributions = groups
            .into_iter()
            .collect::<Vec<_>>()
            .into_par_iter()
            .filter_map(|((flat_index, coefficient_class, left_negative), atoms)| {
                let mut contribution = match atoms.len() {
                    0 => Atom::Zero,
                    1 => atoms
                        .into_iter()
                        .next()
                        .expect("single grouped atom exists"),
                    _ => Atom::add_many(&atoms),
                };
                if left_negative {
                    contribution = -contribution;
                }
                contribution = multiply_atom_by_numeric_coefficient(
                    contribution,
                    &self.coefficient_classes[coefficient_class],
                );
                (!contribution.as_view().is_zero()).then_some((flat_index, contribution))
            })
            .collect::<Vec<_>>();

        Ok(ParamTensor::composite(DataTensor::Sparse(SparseTensor {
            elements: self.collect_output_groups(contributions),
            zero: Atom::Zero,
            structure: self.result_structure.clone(),
        })))
    }

    fn collect_output_groups(
        &self,
        contributions: Vec<(FlatIndex, Atom)>,
    ) -> HashMap<FlatIndex, Atom> {
        let mut output_groups = HashMap::<FlatIndex, Vec<Atom>>::new();
        for (flat_index, contribution) in contributions {
            output_groups
                .entry(flat_index)
                .or_default()
                .push(contribution);
        }

        output_groups
            .into_iter()
            .collect::<Vec<_>>()
            .into_par_iter()
            .filter_map(|(flat_index, atoms)| {
                let atom = match atoms.len() {
                    0 => Atom::Zero,
                    1 => atoms.into_iter().next().expect("single output atom exists"),
                    _ => Atom::add_many(&atoms),
                };
                (!atom.as_view().is_zero()).then_some((flat_index, atom))
            })
            .collect::<HashMap<_, _>>()
    }
}

#[cfg(feature = "shadowing")]
fn common_sparse_atom_factors<S>(sum_sparse: &SparseTensor<Atom, S>) -> Vec<Atom>
where
    S: TensorStructure,
{
    let mut entries = sum_sparse
        .iter_expanded()
        .filter(|(_, atom)| !atom.as_view().is_zero());
    let Some((_, first)) = entries.next() else {
        return Vec::new();
    };

    let mut common = atom_factors(first.as_view())
        .into_iter()
        .filter(|factor| !matches!(factor.as_view(), AtomView::Num(_)) && factor != &Atom::num(1))
        .collect::<Vec<_>>();

    for (_, atom) in entries {
        common.retain(|factor| atom_contains_factor(atom.as_view(), factor));
        if common.is_empty() {
            break;
        }
    }

    common
}

#[cfg(feature = "shadowing")]
fn atom_factors(atom: AtomView<'_>) -> Vec<Atom> {
    match atom {
        AtomView::Mul(mul) => mul
            .iter()
            .map(|factor| factor.as_atom_view().to_owned())
            .collect(),
        _ => vec![atom.to_owned()],
    }
}

#[cfg(feature = "shadowing")]
fn atom_contains_factor(atom: AtomView<'_>, factor: &Atom) -> bool {
    match atom {
        AtomView::Mul(mul) => mul
            .iter()
            .any(|candidate| candidate.as_atom_view() == factor.as_view()),
        _ => atom == factor.as_view(),
    }
}

#[cfg(feature = "shadowing")]
fn atom_without_factors(atom: AtomView<'_>, factors: &[Atom]) -> Atom {
    match atom {
        AtomView::Mul(mul) => {
            let mut used = vec![false; factors.len()];
            let kept = mul
                .iter()
                .filter_map(|candidate| {
                    let candidate = candidate.as_atom_view();
                    if let Some((position, _)) =
                        factors.iter().enumerate().find(|(position, factor)| {
                            !used[*position] && candidate == factor.as_view()
                        })
                    {
                        used[position] = true;
                        None
                    } else {
                        Some(candidate.to_owned())
                    }
                })
                .collect::<Vec<_>>();
            atom_product(kept)
        }
        _ if factors.len() == 1 && atom == factors[0].as_view() => Atom::num(1),
        _ => atom.to_owned(),
    }
}

#[cfg(feature = "shadowing")]
fn atom_product(factors: impl IntoIterator<Item = Atom>) -> Atom {
    factors
        .into_iter()
        .fold(Atom::num(1), |product, factor| product * factor)
}

#[cfg(feature = "shadowing")]
fn collect_numeric_contract_entries<S>(
    numeric: &SparseTensor<Atom, S>,
    numeric_match_axes: &[usize],
    left_match_axes: Option<&[usize]>,
) -> Result<Vec<NumericContractEntry>, ContractionError>
where
    S: TensorStructure,
{
    let mut entries = Vec::with_capacity(numeric.elements.len());
    for (index, coefficient) in numeric.iter_expanded() {
        if coefficient.as_view().is_zero() {
            continue;
        }
        let left_metric_negative = if let Some(left_match_axes) = left_match_axes {
            metric_negative(&numeric.structure, left_match_axes, &index.indices)?
        } else {
            false
        };
        if !numeric_match_axes
            .iter()
            .all(|axis| index.indices.get(*axis).is_some())
        {
            return Err(eyre!("numeric contraction axis out of bounds").into());
        }
        entries.push(NumericContractEntry {
            index: index.indices,
            coefficient: coefficient.clone(),
            coefficient_class: 0,
            left_metric_negative,
        });
    }
    Ok(entries)
}

#[cfg(feature = "shadowing")]
fn component_key(index: &[ConcreteIndex], axes: &[usize]) -> Vec<ConcreteIndex> {
    axes.iter().map(|axis| index[*axis]).collect()
}

#[cfg(feature = "shadowing")]
fn match_flags(order: usize, axes: &[usize]) -> Vec<bool> {
    let mut flags = vec![false; order];
    for axis in axes {
        flags[*axis] = true;
    }
    flags
}

#[cfg(feature = "shadowing")]
fn metric_negative<S>(
    structure: &S,
    matched_axes: &[usize],
    index: &[ConcreteIndex],
) -> Result<bool, ContractionError>
where
    S: TensorStructure,
{
    let slots = structure.external_structure();
    let mut negative = false;
    for axis in matched_axes {
        let slot = slots
            .get(*axis)
            .ok_or_else(|| eyre!("matched axis out of bounds"))?;
        let component = *index
            .get(*axis)
            .ok_or_else(|| eyre!("matched component out of bounds"))?;
        negative ^= *slot
            .rep()
            .negative()?
            .get(component)
            .ok_or_else(|| eyre!("matched component out of representation bounds"))?;
    }
    Ok(negative)
}

#[cfg(feature = "shadowing")]
fn contract_result_index<S>(
    result_structure: &S,
    left_slots: &[S::Slot],
    right_slots: &[S::Slot],
    left_match_flags: &[bool],
    right_match_flags: &[bool],
    left_index: &[ConcreteIndex],
    right_index: &[ConcreteIndex],
) -> Result<Vec<ConcreteIndex>, ContractionError>
where
    S: TensorStructure,
{
    result_structure
        .external_structure_iter()
        .map(|slot| {
            if let Some((position, _)) = left_slots
                .iter()
                .enumerate()
                .find(|(position, left_slot)| !left_match_flags[*position] && **left_slot == slot)
            {
                return left_index
                    .get(position)
                    .copied()
                    .ok_or_else(|| eyre!("left result component out of bounds"));
            }
            if let Some((position, _)) =
                right_slots
                    .iter()
                    .enumerate()
                    .find(|(position, right_slot)| {
                        !right_match_flags[*position] && **right_slot == slot
                    })
            {
                return right_index
                    .get(position)
                    .copied()
                    .ok_or_else(|| eyre!("right result component out of bounds"));
            }
            Err(eyre!(
                "result slot does not map to either contraction operand"
            ))
        })
        .collect::<Result<Vec<_>, _>>()
        .map_err(ContractionError::from)
}

#[cfg(feature = "shadowing")]
fn multiply_atom_by_numeric_coefficient(atom: Atom, coefficient: &Atom) -> Atom {
    if coefficient == &Atom::num(1) {
        atom
    } else if coefficient == &Atom::num(-1) {
        -atom
    } else {
        atom * coefficient.clone()
    }
}

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    PartialEq,
    Eq,
)]
#[cfg_attr(
    feature = "shadowing",
    trait_decode(trait = symbolica::state::HasStateMap),
)]
pub struct Network<S, LibKey, FunKey, Aind = AbstractIndex> {
    pub graph: NetworkGraph<LibKey, FunKey, Aind>,
    pub store: S,
    pub state: NetworkState,
}

#[cfg(feature = "shadowing")]
#[derive(Debug, Clone)]
pub struct ScalarAliases {
    aliased_originals: Vec<bool>,
    aliases_created: usize,
    aliased_terms: usize,
    aliased_bytes: usize,
    max_aliased_bytes: usize,
}

#[cfg(feature = "shadowing")]
impl ScalarAliases {
    fn new(scalars: usize) -> Self {
        Self {
            aliased_originals: vec![false; scalars],
            aliases_created: 0,
            aliased_terms: 0,
            aliased_bytes: 0,
            max_aliased_bytes: 0,
        }
    }

    fn observe_alias(&mut self, index: usize, scalar: &Atom) {
        let view = scalar.as_view();
        let bytes = view.get_byte_size();
        self.aliased_originals[index] = true;
        self.aliases_created += 1;
        self.aliased_terms += scalar.nterms();
        self.aliased_bytes += bytes;
        self.max_aliased_bytes = self.max_aliased_bytes.max(bytes);
    }

    pub fn aliases_created(&self) -> usize {
        self.aliases_created
    }

    pub fn aliased_terms(&self) -> usize {
        self.aliased_terms
    }

    pub fn aliased_bytes(&self) -> usize {
        self.aliased_bytes
    }

    pub fn max_aliased_bytes(&self) -> usize {
        self.max_aliased_bytes
    }

    pub fn is_empty(&self) -> bool {
        self.aliases_created == 0
    }

    pub fn is_aliased(&self, index: usize) -> bool {
        self.aliased_originals.get(index).copied().unwrap_or(false)
    }

    pub fn aliased_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.aliased_originals
            .iter()
            .enumerate()
            .filter_map(|(index, is_aliased)| is_aliased.then_some(index))
    }

    pub fn to_aliased_atom<T>(&self, store: &NetworkStore<T, Atom>, root: Atom) -> AliasedAtom {
        let mut aliased = AliasedAtom::from(root);
        for index in self.aliased_indices() {
            let original = store
                .scalar
                .get(index)
                .expect("scalar alias references an existing scalar")
                .clone();
            aliased = aliased.add_alias(scalar_store_alias(index), original);
        }
        aliased
    }

    pub fn resolve_atom<T>(&self, store: &NetworkStore<T, Atom>, root: Atom) -> Atom {
        if self.is_empty() {
            return root;
        }
        self.to_aliased_atom(store, root).into_inner()
    }
}

#[cfg(feature = "shadowing")]
impl<T, K: Debug, FK: Debug, Aind: AbsInd> Network<NetworkStore<T, Atom>, K, FK, Aind> {
    pub fn alias_scalar_refs(
        &mut self,
        mut should_alias: impl FnMut(usize, &Atom) -> bool,
    ) -> ScalarAliases {
        let mut aliases = ScalarAliases::new(self.store.scalar.len());
        self.store
            .scalar_aliases
            .resize_with(self.store.scalar.len(), || None);

        for (index, scalar) in self.store.scalar.iter().enumerate() {
            if should_alias(index, scalar) {
                aliases.observe_alias(index, scalar);
                self.store.scalar_aliases[index] = Some(scalar_store_alias(index));
            } else {
                self.store.scalar_aliases[index] = None;
            }
        }

        self.graph
            .map_scalar_refs(|scalar| scalar.with_alias_state(|index| aliases.is_aliased(index)));

        aliases
    }

    pub fn resolve_scalar_aliases(&self, aliases: &ScalarAliases, root: Atom) -> Atom {
        aliases.resolve_atom(&self.store, root)
    }

    pub fn aliased_atom(&self, aliases: &ScalarAliases, root: Atom) -> AliasedAtom {
        aliases.to_aliased_atom(&self.store, root)
    }
}

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    PartialEq,
    Eq,
    Copy,
)]
pub enum NetworkState {
    PureScalar,
    Tensor,
    SelfDualTensor,
    Scalar,
}

impl NetworkState {
    pub fn is_scalar(&self) -> bool {
        matches!(self, NetworkState::Scalar | NetworkState::PureScalar)
    }

    pub fn is_tensor(&self) -> bool {
        matches!(self, NetworkState::Tensor | NetworkState::SelfDualTensor)
    }

    pub fn pow(self, pow: i8) -> Self {
        match self {
            NetworkState::PureScalar => NetworkState::PureScalar,
            NetworkState::Scalar => NetworkState::Scalar,
            NetworkState::SelfDualTensor => {
                if pow % 2 == 0 {
                    NetworkState::Scalar
                } else {
                    NetworkState::SelfDualTensor
                }
            }
            NetworkState::Tensor => panic!("Cannot have integer power of non-self dual tensor"),
        }
    }

    pub fn is_compatible(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (NetworkState::Tensor, NetworkState::Tensor)
                | (NetworkState::SelfDualTensor, NetworkState::SelfDualTensor)
                | (NetworkState::Scalar, NetworkState::Scalar)
                | (NetworkState::PureScalar, NetworkState::Scalar)
                | (NetworkState::Scalar, NetworkState::PureScalar)
                | (NetworkState::PureScalar, NetworkState::PureScalar)
        )
    }
}
impl MulAssign for NetworkState {
    fn mul_assign(&mut self, rhs: Self) {
        // println!("{self:?} *={rhs:?}");
        *self = match (*self, rhs) {
            (NetworkState::PureScalar, NetworkState::PureScalar) => NetworkState::PureScalar,
            (NetworkState::PureScalar, NetworkState::Scalar) => NetworkState::Scalar,
            (NetworkState::Tensor, _) => NetworkState::Tensor,
            (NetworkState::Scalar, NetworkState::PureScalar) => NetworkState::Scalar,
            (_, NetworkState::Tensor) => NetworkState::Tensor,
            (NetworkState::SelfDualTensor, _) => NetworkState::SelfDualTensor,
            (_, NetworkState::SelfDualTensor) => NetworkState::SelfDualTensor,
            (NetworkState::Scalar, NetworkState::Scalar) => NetworkState::Scalar,
        }
    }
}

impl AddAssign for NetworkState {
    fn add_assign(&mut self, rhs: Self) {
        // println!("{self:?} *={rhs:?}");
        *self = match (*self, rhs) {
            (NetworkState::PureScalar, NetworkState::PureScalar) => NetworkState::PureScalar,
            (NetworkState::PureScalar, NetworkState::Scalar) => NetworkState::Scalar,
            (a, b) => {
                assert_eq!(a, b, "Cannot add incompatible network states:{a:?} + {b:?}");
                a
            }
        }
    }
}

// pub type TensorNetwork<T, S, Str: TensorScalarStore<Tensor = T, Scalar = S>, K> = Network<Str, K>;

// pub struct TensorNetwork<
//     T,
//     S,
//     K,
//     Str: TensorScalarStore<Tensor = T, Scalar = S> = NetworkStore<T, S>,
// > {
//     net: Network<Str, K>,
// }

pub mod graph;
pub mod library;
#[doc(hidden)]
pub mod profile;
pub mod set;
pub mod store;

impl<S: TensorScalarStoreMapping, K: Clone, FK: Clone, Aind: AbsInd> TensorScalarStoreMapping
    for Network<S, K, FK, Aind>
{
    type Store<U, V> = Network<S::Store<U, V>, K, FK, Aind>;
    type Scalar = S::Scalar;
    type Tensor = S::Tensor;

    fn iter_scalars(&self) -> impl Iterator<Item = &Self::Scalar> {
        self.store.iter_scalars()
    }

    fn iter_tensors(&self) -> impl Iterator<Item = &Self::Tensor> {
        self.store.iter_tensors()
    }

    fn iter_scalars_mut(&mut self) -> impl Iterator<Item = &mut Self::Scalar> {
        self.store.iter_scalars_mut()
    }
    fn iter_tensors_mut(&mut self) -> impl Iterator<Item = &mut Self::Tensor> {
        self.store.iter_tensors_mut()
    }

    fn map<U, V>(
        self,
        scalar_map: impl FnMut(Self::Scalar) -> U,
        tensor_map: impl FnMut(Self::Tensor) -> V,
    ) -> Self::Store<V, U> {
        Network {
            store: self.store.map(scalar_map, tensor_map),
            graph: self.graph,
            state: self.state,
        }
    }

    fn map_result<U, V, Er>(
        self,
        scalar_map: impl FnMut(Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er> {
        Ok(Network {
            store: self.store.map_result(scalar_map, tensor_map)?,
            graph: self.graph,
            state: self.state,
        })
    }

    fn map_ref<'a, U, V>(
        &'a self,
        scalar_map: impl FnMut(&'a Self::Scalar) -> U,
        tensor_map: impl FnMut(&'a Self::Tensor) -> V,
    ) -> Self::Store<V, U> {
        Network {
            store: self.store.map_ref(scalar_map, tensor_map),
            graph: self.graph.clone(),
            state: self.state,
        }
    }

    fn map_ref_result<U, V, Er>(
        &self,
        scalar_map: impl FnMut(&Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(&Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er> {
        Ok(Network {
            store: self.store.map_ref_result(scalar_map, tensor_map)?,
            graph: self.graph.clone(),
            state: self.state,
        })
    }

    fn map_ref_enumerate<U, V>(
        &self,
        scalar_map: impl FnMut((usize, &Self::Scalar)) -> U,
        tensor_map: impl FnMut((usize, &Self::Tensor)) -> V,
    ) -> Self::Store<V, U> {
        Network {
            store: self.store.map_ref_enumerate(scalar_map, tensor_map),
            graph: self.graph.clone(),
            state: self.state,
        }
    }

    fn map_ref_result_enumerate<U, V, Er>(
        &self,
        scalar_map: impl FnMut((usize, &Self::Scalar)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &Self::Tensor)) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er> {
        Ok(Network {
            store: self
                .store
                .map_ref_result_enumerate(scalar_map, tensor_map)?,
            graph: self.graph.clone(),
            state: self.state,
        })
    }

    fn map_ref_mut<U, V>(
        &mut self,
        scalar_map: impl FnMut(&mut Self::Scalar) -> U,
        tensor_map: impl FnMut(&mut Self::Tensor) -> V,
    ) -> Self::Store<V, U> {
        Network {
            store: self.store.map_ref_mut(scalar_map, tensor_map),
            graph: self.graph.clone(),
            state: self.state,
        }
    }

    fn map_ref_mut_result<U, V, Er>(
        &mut self,
        scalar_map: impl FnMut(&mut Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(&mut Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er> {
        Ok(Network {
            store: self.store.map_ref_mut_result(scalar_map, tensor_map)?,
            graph: self.graph.clone(),
            state: self.state,
        })
    }

    fn map_ref_mut_enumerate<U, V>(
        &mut self,
        scalar_map: impl FnMut((usize, &mut Self::Scalar)) -> U,
        tensor_map: impl FnMut((usize, &mut Self::Tensor)) -> V,
    ) -> Self::Store<V, U> {
        Network {
            store: self.store.map_ref_mut_enumerate(scalar_map, tensor_map),
            graph: self.graph.clone(),
            state: self.state,
        }
    }

    fn map_ref_mut_result_enumerate<U, V, Er>(
        &mut self,
        scalar_map: impl FnMut((usize, &mut Self::Scalar)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &mut Self::Tensor)) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er> {
        Ok(Network {
            store: self
                .store
                .map_ref_mut_result_enumerate(scalar_map, tensor_map)?,
            graph: self.graph.clone(),
            state: self.state,
        })
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> Default for Network<S, K, FK, Aind> {
    fn default() -> Self {
        Self::one()
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> NMul for Network<S, K, FK, Aind> {
    type Output = Self;
    fn n_mul<I: IntoIterator<Item = Self>>(self, iter: I) -> Self::Output {
        let _span = profile::span(Timer::NetworkNMul);
        profile::bump(Counter::NetworkNMul, 1);
        let mut store = self.store;
        let mut state = self.state;

        let items = iter.into_iter().map(|mut a| {
            let _span = profile::span(Timer::NetworkNMulStorePrep);
            a.graph.shift_scalars(store.n_scalars());
            a.graph.shift_tensors(store.n_tensors());
            store.extend(a.store);

            state *= a.state;
            a.graph
        });

        let graph = {
            let _span = profile::span(Timer::NetworkNMulGraph);
            self.graph.n_mul(items)
        };

        if state.is_tensor() {
            let _span = profile::span(Timer::NetworkNMulDangling);
            state = graph.state();
        }

        Network {
            graph,
            store,
            state,
        }
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> Mul for Network<S, K, FK, Aind> {
    type Output = Self;
    fn mul(mut self, mut other: Self) -> Self::Output {
        let mut store = self.store;

        other.graph.shift_scalars(store.n_scalars());
        other.graph.shift_tensors(store.n_tensors());
        store.extend(other.store);
        self.state *= other.state;

        Network {
            graph: self.graph * other.graph,
            store,
            state: self.state,
        }
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> MulAssign
    for Network<S, K, FK, Aind>
{
    fn mul_assign(&mut self, mut rhs: Self) {
        rhs.graph.shift_scalars(self.store.n_scalars());
        rhs.graph.shift_tensors(self.store.n_tensors());
        self.store.extend(rhs.store);
        self.state *= rhs.state;
        self.graph *= rhs.graph;
    }
}

impl<T: TensorStructure, S, FK: Debug, K: Debug, Aind: AbsInd> MulAssign<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    fn mul_assign(&mut self, rhs: T) {
        *self *= Network::from_tensor(rhs);
    }
}

impl<T: TensorStructure, S, FK: Debug, K: Debug, Aind: AbsInd> Mul<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    type Output = Self;
    fn mul(self, other: T) -> Self::Output {
        let mut store = self.store;

        let mut other = Network::from_tensor(other);

        other.graph.shift_scalars(store.n_scalars());
        other.graph.shift_tensors(store.n_tensors());
        store.extend(other.store);

        Network {
            graph: self.graph * other.graph,
            store,
            state: self.state,
        }
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> Add for Network<S, K, FK, Aind> {
    type Output = Self;
    fn add(self, mut other: Self) -> Self::Output {
        let mut store = self.store;

        other.graph.shift_scalars(store.n_scalars());
        other.graph.shift_tensors(store.n_tensors());
        store.extend(other.store);

        Network {
            graph: self.graph + other.graph,
            store,
            state: self.state,
        }
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> AddAssign
    for Network<S, K, FK, Aind>
{
    fn add_assign(&mut self, mut rhs: Self) {
        rhs.graph.shift_scalars(self.store.n_scalars());
        rhs.graph.shift_tensors(self.store.n_tensors());
        self.store.extend(rhs.store);

        self.graph += rhs.graph;
    }
}

impl<T: TensorStructure, S, FK: Debug, K: Debug, Aind: AbsInd> AddAssign<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    fn add_assign(&mut self, rhs: T) {
        *self += Network::from_tensor(rhs);
    }
}

impl<T: TensorStructure, S, FK: Debug, K: Debug, Aind: AbsInd> Add<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    type Output = Self;
    fn add(mut self, other: T) -> Self::Output {
        self += other;
        self
    }
}

impl<T: TensorStructure, FK: Debug, K: Debug, Aind: AbsInd> Add<i8>
    for Network<NetworkStore<T, i8>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    type Output = Self;
    fn add(mut self, other: i8) -> Self::Output {
        let mut other = Network::from_scalar(other);
        other.graph.shift_tensors(self.store.n_tensors());
        other.graph.shift_tensors(self.store.n_tensors());

        self.store.extend(other.store);
        Network {
            graph: self.graph + other.graph,
            store: self.store,
            state: self.state,
        }
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> NAdd for Network<S, K, FK, Aind> {
    type Output = Self;
    fn n_add<I: IntoIterator<Item = Self>>(self, iter: I) -> Self::Output {
        let _span = profile::span(Timer::NetworkNAdd);
        profile::bump(Counter::NetworkNAdd, 1);
        let mut store = self.store;

        let mut state = self.state;

        let items = iter.into_iter().map(|mut a| {
            let _span = profile::span(Timer::NetworkNAddStorePrep);
            a.graph.shift_scalars(store.n_scalars());
            a.graph.shift_tensors(store.n_tensors());
            store.extend(a.store);

            state += a.state;
            a.graph
        });

        Network {
            graph: {
                let _span = profile::span(Timer::NetworkNAddGraph);
                self.graph.n_add(items)
            },
            store,
            state,
        }
    }
}

impl<S: TensorScalarStore, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Neg
    for Network<S, K, FK, Aind>
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            store: self.store,
            graph: self.graph.neg(),
            state: self.state,
        }
    }
}

impl<S: TensorScalarStore, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Sub
    for Network<S, K, FK, Aind>
{
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<S: TensorScalarStore, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> SubAssign
    for Network<S, K, FK, Aind>
{
    fn sub_assign(&mut self, mut rhs: Self) {
        rhs.graph.shift_scalars(self.store.n_scalars());
        rhs.graph.shift_tensors(self.store.n_tensors());
        self.store.extend(rhs.store);

        self.graph -= rhs.graph
    }
}

impl<T: TensorStructure, S, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> SubAssign<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    fn sub_assign(&mut self, rhs: T) {
        *self -= Network::from_tensor(rhs)
    }
}

impl<T: TensorStructure, S, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Sub<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    type Output = Self;
    fn sub(mut self, other: T) -> Self::Output {
        self -= other;
        self
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> Network<S, K, FK, Aind> {
    pub fn pow(self, pow: i8) -> Self {
        Self {
            store: self.store,
            graph: self.graph.pow(pow),
            state: self.state.pow(pow),
        }
    }
    pub fn fun(self, key: FK) -> Self {
        let graph = self.graph.function(key);
        Self {
            store: self.store,
            state: graph.state(),
            graph,
        }
    }

    pub fn from_scalar(scalar: S::Scalar) -> Self {
        let _span = profile::span(Timer::FromScalar);
        profile::bump(Counter::FromScalar, 1);
        let mut store = S::default();
        let id = store.add_scalar(scalar);
        Network {
            graph: NetworkGraph::scalar(id),
            store,

            state: NetworkState::PureScalar,
        }
    }

    pub fn merge_ops(&mut self)
    where
        K: Clone,
    {
        self.graph.merge_ops();
    }

    pub fn from_tensor(tensor: S::Tensor) -> Self
    where
        S::Tensor: TensorStructure,
        <S::Tensor as TensorStructure>::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let _span = profile::span(Timer::FromTensor);
        profile::bump(Counter::FromTensor, 1);
        let mut store = S::default();

        let id = store.add_tensor(tensor);
        let graph = NetworkGraph::tensor(store.get_tensor(id), NetworkLeaf::LocalTensor(id));
        let state = graph.state();

        Network {
            graph,
            store,
            state,
        }
    }

    pub fn library_tensor<T>(tensor: &T, key: PermutedStructure<K>) -> Self
    where
        T: TensorStructure,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let _span = profile::span(Timer::LibraryTensor);
        profile::bump(Counter::LibraryTensor, 1);
        let indices = tensor.external_indices_iter().collect();
        let graph = NetworkGraph::tensor(tensor, NetworkLeaf::LibraryKey { key, indices });
        let state = graph.state();
        Network {
            graph,
            store: S::default(),
            state,
        }
    }

    pub fn one() -> Self {
        Network {
            graph: NetworkGraph::one(),
            store: S::default(),
            state: NetworkState::PureScalar,
        }
    }

    pub fn zero() -> Self {
        Network {
            graph: NetworkGraph::zero(),
            store: S::default(),
            state: NetworkState::PureScalar,
        }
    }
}

#[derive(Error, Debug)]
pub enum TensorNetworkError<K: Display, FK: Display> {
    #[error("Slot edge to prod node")]
    SlotEdgeToProdNode,
    #[error("Slot edge to scalar node")]
    SlotEdgeToScalarNode,
    #[error("More than one neg")]
    MoreThanOneNeg,
    #[error("Childless neg")]
    ChildlessNeg,
    #[error("Contraction Error:{0}")]
    ContractionError(#[from] ContractionError),
    #[error("Scalar connected by a slot edge")]
    ScalarSlotEdge,
    #[error("Structure Error:{0}")]
    StructErr(#[from] StructureError),
    #[error("LibraryError:{0}")]
    LibErr(#[from] LibraryError<K>),
    #[error("FunctionLibraryError:{0}")]
    FunLibErr(#[from] FunctionLibraryError<FK>),
    #[error("Non tensor node still present")]
    NonTensorNodePresent,
    #[error("Negative non-even power on non-scalar node:{0}")]
    NegativeExponentNonScalar(String),
    #[error("Too many arguments for function:{0}")]
    TooManyArgsFunction(String),
    #[error("Non self-dual tensor power{0}")]
    InvalidDotFunction(String),
    #[error("Invalid dot function{0}")]
    NonSelfDualTensorPower(String),
    #[error("invalid resulting node{0}")]
    InvalidResultNode(NetworkNode<DummyKey, FK>),
    #[error("internal edge still present, contract it first")]
    InternalEdgePresent,
    #[error("uncontracted scalar")]
    UncontractedScalar,
    #[error("Cannot contract edge between {0} and {1}")]
    CannotContractEdgeBetween(String, String),
    #[error("no nodes in the graph")]
    NoNodes,
    #[error("no scalar present")]
    NoScalar,
    #[error("more than one node in the graph")]
    MoreThanOneNode,
    #[error("is not scalar output")]
    NotScalarOutput,
    #[error("failed scalar multiplication")]
    FailedScalarMul,
    #[error("scalar field is empty")]
    ScalarFieldEmpty,
    #[error("not all scalars: {0}")]
    NotAllScalars(String),
    #[error("try to sum scalar with library tensor: {0}")]
    ScalarLibSum(String),
    #[error("try to sum scalar with a tensor: {0}")]
    SumScalarTensor(String),
    #[error("Incompatible summands: {0}")]
    IncompatibleSummand(String),
    #[error("failed to contract")]
    FailedContract(ContractionError),
    #[error("negative exponent not yet supported")]
    NegativeExponent,
    #[error("failed to contract: {0}")]
    FailedContractMsg(String),
    #[error(transparent)]
    Other(#[from] eyre::Error),
    #[error("Io error")]
    InOut(#[from] std::io::Error),
    #[error("Infallible")]
    Infallible,
}

impl<K: Display, FK: Display> From<Infallible> for TensorNetworkError<K, FK> {
    fn from(_: Infallible) -> Self {
        TensorNetworkError::Infallible
    }
}

pub enum TensorOrScalarOrKey<T, S, K, Aind> {
    Tensor {
        tensor: T,
        graph_slots: Vec<LibrarySlot<Aind>>,
    },
    Scalar(S),
    Key {
        key: K,
        nodeid: NodeIndex,
    },
}

pub enum ExecutionResult<T> {
    One,
    Zero,
    Val(T),
}

impl<T: Display> Display for ExecutionResult<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ExecutionResult::One => write!(f, "One"),
            ExecutionResult::Zero => write!(f, "Zero"),
            ExecutionResult::Val(val) => write!(f, "{}", val),
        }
    }
}

impl<
    T: TensorStructure,
    S,
    K: Display + Debug,
    FK: Display + Debug,
    Str: TensorScalarStore<Tensor = T, Scalar = S>,
    Aind: AbsInd,
> Network<Str, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    pub fn validate(&self)
    where
        K: TensorStructure,
    {
        for (n, _neigh, v) in self.graph.graph.iter_nodes() {
            match v {
                NetworkNode::Leaf(NetworkLeaf::LibraryKey { key, .. }) => {
                    let reps = self
                        .graph
                        .slots(n)
                        .into_iter()
                        .map(|s| s.rep())
                        .collect::<Vec<_>>();
                    // let p = Permutation::sort(&reps);

                    let n_reps = key
                        .structure
                        .external_reps_iter()
                        .map(|r| r.to_lib())
                        .collect::<Vec<_>>();
                    // let q = Permutation::sort(&n_reps);
                    // println!("p{p}q{q}");
                    assert_eq!(n_reps, reps);
                }
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(k)) => {
                    let reps = self
                        .graph
                        .slots(n)
                        .into_iter()
                        .map(|s| s.rep())
                        .collect::<Vec<_>>();
                    let n_reps = self
                        .store
                        .get_tensor(*k)
                        .external_reps_iter()
                        .map(|r| r.to_lib())
                        .collect::<Vec<_>>();
                    assert_eq!(n_reps, reps);
                }
                _ => {}
            }
        }
    }

    #[allow(clippy::result_large_err, clippy::type_complexity)]
    pub fn result(
        &self,
    ) -> Result<
        ExecutionResult<TensorOrScalarOrKey<&T, &S, &PermutedStructure<K>, Aind>>,
        TensorNetworkError<K, FK>,
    >
    where
        FK: Clone,
    {
        let (node, nid, graph_slots) = self.graph.result()?;

        match node {
            NetworkNode::Leaf(l) => match l {
                NetworkLeaf::LibraryKey { key, .. } => {
                    Ok(ExecutionResult::Val(TensorOrScalarOrKey::Key {
                        key,
                        nodeid: nid,
                    }))
                }
                NetworkLeaf::LocalTensor(t) => {
                    Ok(ExecutionResult::Val(TensorOrScalarOrKey::Tensor {
                        tensor: self.store.get_tensor(*t),
                        graph_slots,
                    }))
                }
                NetworkLeaf::TensorSum(_) => Err(TensorNetworkError::Other(eyre!(
                    "result is an unmaterialized tensor sum"
                ))),
                NetworkLeaf::TensorTerm(_) => Err(TensorNetworkError::Other(eyre!(
                    "result is an unmaterialized scaled tensor"
                ))),
                NetworkLeaf::TensorTermSum(_) => Err(TensorNetworkError::Other(eyre!(
                    "result is an unmaterialized scaled tensor sum"
                ))),
                NetworkLeaf::Scalar(t) => Ok(ExecutionResult::Val(TensorOrScalarOrKey::Scalar(
                    self.store.get_scalar_ref(*t),
                ))),
            },
            NetworkNode::Op(o) => match o {
                NetworkOp::Product => Ok(ExecutionResult::One),
                NetworkOp::Sum => Ok(ExecutionResult::Zero),
                o => Err(TensorNetworkError::InvalidResultNode(NetworkNode::Op(
                    o.clone(),
                ))),
            },
        }
    }

    #[allow(clippy::result_large_err)]
    pub fn result_tensor<'a, LT, L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>>(
        &'a self,
        lib: &L,
    ) -> Result<ExecutionResult<Cow<'a, T>>, TensorNetworkError<K, FK>>
    where
        S: 'a,
        T: Clone
            + ScalarTensor
            + HasStructure
            + Ref
            + FastTensorSum
            + ScalarMul<S, Output = T>
            + for<'b> AddAssign<<T as Ref>::Ref<'b>>,
        S: Clone + Into<T::Scalar>,
        K: Display + Debug,
        FK: Display + Debug + Clone,
        LT: TensorStructure<Indexed = T> + Clone + LibraryTensor<WithIndices = T>,
        T: PermuteTensor<Permuted = T>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        let tensor_term_owned = |term: &TensorTerm| -> Result<T, TensorNetworkError<K, FK>> {
            let tensor = self.store.get_tensor(term.tensor);
            match term.scalar {
                Some(scalar) => tensor
                    .scalar_mul(self.store.get_scalar_ref(scalar))
                    .ok_or(TensorNetworkError::FailedScalarMul),
                None => Ok(tensor.clone()),
            }
        };

        let tensor_terms_owned = |terms: &[TensorTerm]| -> Result<T, TensorNetworkError<K, FK>> {
            let materialized = terms
                .iter()
                .map(tensor_term_owned)
                .collect::<Result<Vec<_>, _>>()?;
            if let Some(result) = {
                let refs = materialized.iter().collect::<Vec<_>>();
                T::fast_tensor_sum(&refs, None)
            } {
                return Ok(result);
            }
            Ok(balanced_ref_sum(materialized, None))
        };

        let (node, nodeid, _) = self.graph.result()?;
        if let NetworkNode::Leaf(leaf) = node {
            match leaf {
                NetworkLeaf::TensorSum(indices) => {
                    let terms = indices
                        .iter()
                        .copied()
                        .map(TensorTerm::tensor)
                        .collect::<Vec<_>>();
                    return Ok(ExecutionResult::Val(Cow::Owned(tensor_terms_owned(
                        &terms,
                    )?)));
                }
                NetworkLeaf::TensorTerm(term) => {
                    return Ok(ExecutionResult::Val(Cow::Owned(tensor_term_owned(term)?)));
                }
                NetworkLeaf::TensorTermSum(terms) => {
                    return Ok(ExecutionResult::Val(Cow::Owned(tensor_terms_owned(terms)?)));
                }
                NetworkLeaf::LibraryKey { .. } => {
                    let less = self.graph.get_lib_data(lib, nodeid).unwrap();
                    return Ok(ExecutionResult::Val(Cow::Owned(less)));
                }
                NetworkLeaf::LocalTensor(_) | NetworkLeaf::Scalar(_) => {}
            }
        }

        Ok(match self.result()? {
            ExecutionResult::One => ExecutionResult::One,
            ExecutionResult::Zero => ExecutionResult::Zero,
            ExecutionResult::Val(v) => ExecutionResult::Val(match v {
                TensorOrScalarOrKey::Tensor { tensor, .. } => Cow::Borrowed(tensor),
                TensorOrScalarOrKey::Scalar(s) => Cow::Owned(T::new_scalar(s.clone().into())),
                TensorOrScalarOrKey::Key { nodeid, .. } => {
                    let less = self.graph.get_lib_data(lib, nodeid).unwrap();

                    Cow::Owned(less)
                }
            }),
        })
    }

    #[allow(clippy::result_large_err)]
    pub fn result_scalar<'a>(
        &'a self,
    ) -> Result<ExecutionResult<Cow<'a, S>>, TensorNetworkError<K, FK>>
    where
        T: Clone + ScalarTensor + 'a,
        T::Scalar: Into<S>,
        K: Display,
        FK: Display + Clone,
        S: Clone
            + Ref
            + for<'b> AddAssign<<S as Ref>::Ref<'b>>
            + for<'b> MulAssign<<S as Ref>::Ref<'b>>,
    {
        let tensor_scalar = |tensor: usize| -> Result<S, TensorNetworkError<K, FK>> {
            self.store
                .get_tensor(tensor)
                .clone()
                .scalar()
                .ok_or(TensorNetworkError::NoScalar)
                .map(Into::into)
        };

        let tensor_term_scalar = |term: &TensorTerm| -> Result<S, TensorNetworkError<K, FK>> {
            let mut scalar: S = tensor_scalar(term.tensor)?;
            if let Some(factor) = term.scalar {
                scalar *= self.store.get_scalar_ref(factor).refer();
            }
            Ok(scalar)
        };

        let tensor_sum_scalar =
            |indices: &[usize]| -> Result<ExecutionResult<Cow<'a, S>>, TensorNetworkError<K, FK>> {
                let mut iter = indices.iter();
                let Some(first) = iter.next() else {
                    return Ok(ExecutionResult::Zero);
                };
                let mut accumulator = tensor_scalar(*first)?;
                for tensor in iter {
                    let term_scalar = tensor_scalar(*tensor)?;
                    accumulator += term_scalar.refer();
                }
                Ok(ExecutionResult::Val(Cow::Owned(accumulator)))
            };

        let tensor_term_sum_scalar = |terms: &[TensorTerm]| -> Result<
            ExecutionResult<Cow<'a, S>>,
            TensorNetworkError<K, FK>,
        > {
            let mut iter = terms.iter();
            let Some(first) = iter.next() else {
                return Ok(ExecutionResult::Zero);
            };
            let mut accumulator = tensor_term_scalar(first)?;
            for term in iter {
                let term_scalar = tensor_term_scalar(term)?;
                accumulator += term_scalar.refer();
            }
            Ok(ExecutionResult::Val(Cow::Owned(accumulator)))
        };

        let (node, _, _) = self.graph.result()?;
        if let NetworkNode::Leaf(leaf) = node {
            match leaf {
                NetworkLeaf::TensorSum(indices) => return tensor_sum_scalar(indices),
                NetworkLeaf::TensorTerm(term) => {
                    return Ok(ExecutionResult::Val(Cow::Owned(tensor_term_scalar(term)?)));
                }
                NetworkLeaf::TensorTermSum(terms) => return tensor_term_sum_scalar(terms),
                _ => {}
            }
        }

        Ok(match self.result()? {
            ExecutionResult::One => ExecutionResult::One,
            ExecutionResult::Zero => ExecutionResult::Zero,
            ExecutionResult::Val(v) => ExecutionResult::Val(match v {
                TensorOrScalarOrKey::Tensor { tensor: t, .. } => Cow::Owned(
                    t.clone()
                        .scalar()
                        .ok_or(TensorNetworkError::NoScalar)?
                        .into(),
                ),
                TensorOrScalarOrKey::Scalar(s) => Cow::Borrowed(s),
                TensorOrScalarOrKey::Key { .. } => return Err(TensorNetworkError::NoScalar),
            }),
        })
    }

    pub fn cast<U>(self) -> Network<Str::Store<U, S>, K, FK, Aind>
    where
        K: Clone,
        FK: Clone,
        T: CastStructure<U> + HasStructure,
        T::Structure: TensorStructure,
        U: HasStructure,
        U::Structure: From<T::Structure> + TensorStructure<Slot = T::Slot>,
    {
        self.map(|a| a, |t| t.cast_structure())
    }
}

pub trait StructureLessDisplay {
    fn display(&self) -> String {
        String::new()
    }
}

impl<S: HasName> StructureLessDisplay for S
where
    S::Args: StructureLessDisplay,
    S::Name: Display,
{
    fn display(&self) -> String {
        format!(
            "{}({})",
            self.name().map(|t| t.to_string()).unwrap_or_default(),
            self.args().map(|t| t.display()).unwrap_or_default()
        )
    }
}

impl<S, K: Display + Debug, FK: Display + Debug, Aind: AbsInd> Network<S, K, FK, Aind> {
    pub fn dot(&self) -> std::string::String {
        self.graph.dot()
    }

    pub fn dot_pretty(&self) -> std::string::String
    where
        S: TensorScalarStore,
        S::Scalar: Display,
        K: StructureLessDisplay,
        S::Tensor: StructureLessDisplay,
    {
        self.graph.dot_impl(
            |i| {
                let ss = &self.store.get_scalar(i);
                format!("{}:{}", i, ss)
            },
            |k| k.display(),
            |t| {
                let tt = &self.store.get_tensor(t);
                tt.display()
            },
            |fk| fk.to_string(),
        )
    }
}

impl<T, S, FK: Debug, K: Debug, Aind: AbsInd> Network<NetworkStore<T, S>, K, FK, Aind> {
    pub fn dot_display_impl(
        &self,
        scalar_disp: impl Fn(&S) -> String,
        library_disp: impl Fn(&K) -> Option<String>,
        tensor_disp: impl Fn(&T) -> String,
        function_disp: impl Fn(&FK) -> String,
    ) -> std::string::String {
        self.graph.graph.dot_impl(
            &self.graph.graph.full_filter(),
            "",
            &|_| None,
            &|e| {
                if let NetworkEdge::Slot(s) = e {
                    Some(format!("label=\"{s}\""))
                } else {
                    None
                }
            },
            &|n| match n {
                NetworkNode::Leaf(l) => match l {
                    NetworkLeaf::LibraryKey { key, .. } => {
                        // if let Ok(v) = lib.get(l) {
                        Some(format!("label = \"L:{}\"", library_disp(&key.structure)?))
                        // } else {
                        // None
                        // }
                    }
                    NetworkLeaf::LocalTensor(l) => Some(format!(
                        "label = \"T:{}\"",
                        tensor_disp(self.store.get_tensor(*l))
                    )),
                    NetworkLeaf::TensorSum(terms) => {
                        Some(format!("label = \"TS:{}\"", terms.len()))
                    }
                    NetworkLeaf::TensorTerm(term) => Some(match term.scalar {
                        Some(scalar) => format!(
                            "label = \"TT:{}*{}\"",
                            tensor_disp(self.store.get_tensor(term.tensor)),
                            scalar.display_with(|index| scalar_disp(self.store.get_scalar(index)))
                        ),
                        None => format!(
                            "label = \"TT:{}\"",
                            tensor_disp(self.store.get_tensor(term.tensor))
                        ),
                    }),
                    NetworkLeaf::TensorTermSum(terms) => {
                        Some(format!("label = \"TTS:{}\"", terms.len()))
                    }
                    NetworkLeaf::Scalar(s) => Some(format!(
                        "label = \"S:{}\"",
                        s.display_with(|index| scalar_disp(self.store.get_scalar(index)))
                    )),
                },
                NetworkNode::Op(o) => {
                    Some(format!("label = \"{}\"", o.display_with(&function_disp)))
                }
            },
        )
        // self.graph.dot()
    }
}

// use log::trace;
#[cfg(feature = "shadowing")]
pub mod parsing;
// use log::trace;
pub mod contract;
pub use contract::{
    ContractScalars, ContractionStrategy, DEFAULT_EXACT_JOIN_LIMIT, MinResultRank,
    MinResultRankWith, PAIR_SCORE_ATOM_AWARE, PAIR_SCORE_ENTRY_AWARE, PAIR_SCORE_RESULT_RANK_ONLY,
    PAIR_SCORE_SPARSE_ATOM_AWARE, ProductContraction, SingleSmallestDegree, SmallestDegree,
    SmallestDegreeIter,
};
pub trait ExecutionStrategy<E, FL, L, K, FK, Aind>
where
    E: ExecuteOp<FL, L, K, FK, Aind>,
{
    /// Run the entire contraction to one leaf.
    #[allow(clippy::result_large_err)]
    fn execute_all<C: ContractionStrategy<E, L, K, FK, Aind>>(
        executor: &mut E,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        fnlib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display;
}

fn first_executable_operation<K, FK, Aind>(
    graph: &mut NetworkGraph<K, FK, Aind>,
) -> Option<NetworkOperation<FK>>
where
    K: Debug,
    FK: Clone + Debug,
    Aind: AbsInd,
{
    graph.cache_expr_tree_roots();
    graph.ready_operation_ref().map(|op_ref| (&op_ref).into())
}

#[allow(clippy::result_large_err)]
fn collapse_operation_subgraph<K, FK, Aind>(
    graph: &mut NetworkGraph<K, FK, Aind>,
    operation: &NetworkOperation<FK>,
    replacement: NetworkLeaf<K, Aind>,
) -> Result<(), TensorNetworkError<K, FK>>
where
    K: Debug + Display,
    FK: Debug + Display,
    Aind: AbsInd,
{
    let mut ignored: SuBitGraph = graph.graph.empty_subgraph();
    graph
        .identify_subgraph_nodes_without_deleting_self_edges(
            operation.subgraph(),
            NetworkNode::Leaf(replacement),
            &mut ignored,
        )
        .ok_or_else(|| {
            TensorNetworkError::Other(eyre!("ready operation subgraph did not contain any nodes"))
        })?;
    graph.finish_deferred_node_identifications();
    if !ignored.is_empty() {
        graph.delete(&ignored);
    }
    Ok(())
}

#[allow(clippy::result_large_err)]
fn execute_operation_in_place<C, E, FL, L, K, FK, Aind>(
    executor: &mut E,
    graph: &mut NetworkGraph<K, FK, Aind>,
    operation: &NetworkOperation<FK>,
    lib: &L,
    fnlib: &FL,
    ignored: &mut SuBitGraph,
) -> Result<bool, TensorNetworkError<K, FK>>
where
    C: ContractionStrategy<E, L, K, FK, Aind>,
    E: ExecuteOp<FL, L, K, FK, Aind>,
    K: Clone + Debug + Display,
    FK: Clone + Debug + Display,
    Aind: AbsInd,
{
    if matches!(operation.op(), NetworkOp::Product)
        && <C as ContractionStrategy<E, L, K, FK, Aind>>::SUPPORTS_PARTIAL_GRAPH_REWRITE
    {
        return C::contract_product_in_place(executor, graph, operation, lib, ignored);
    }

    let replacement = executor.execute::<C>(graph, operation, lib, fnlib)?;
    graph
        .identify_subgraph_nodes_without_deleting_self_edges(
            operation.subgraph(),
            NetworkNode::Leaf(replacement),
            ignored,
        )
        .ok_or_else(|| {
            TensorNetworkError::Other(eyre!("ready operation subgraph did not contain any nodes"))
        })?;
    graph.finish_deferred_node_identifications();
    Ok(true)
}

fn plan_ready_operation_batch<K, FK, Aind>(
    graph: &mut NetworkGraph<K, FK, Aind>,
    ignored: &SuBitGraph,
) -> Vec<NetworkOperation<FK>>
where
    K: Debug,
    FK: Clone + Debug,
    Aind: AbsInd,
{
    let _span = profile::span(Timer::ExecuteFindReady);
    let rerooted = {
        let _span = profile::span(Timer::ExecuteCacheRoots);
        graph.cache_expr_tree_roots_ignoring(ignored)
    };
    let planned = {
        let _span = profile::span(Timer::ExecuteReadyBatch);
        graph.ready_operations_from_tree_ignoring(ignored)
    };
    let batch_len = planned.len();
    let batch_subgraph_hedges = if profile::enabled() {
        planned
            .iter()
            .map(|op| op.subgraph().n_included())
            .sum::<usize>()
    } else {
        0
    };

    if profile::enabled() {
        eprintln!(
            "spenso_profile execute.plan graph_nodes={} graph_hedges={} ignored_hedges={} rerooted={} ready_batch={} batch_subgraph_hedges={}",
            graph.graph.n_nodes(),
            graph.graph.n_hedges(),
            ignored.n_included(),
            rerooted,
            batch_len,
            batch_subgraph_hedges,
        );
    }

    planned
}

fn log_sum_add(
    sum_start: Option<std::time::Instant>,
    done: usize,
    total: usize,
    add_start: std::time::Instant,
) {
    if let Some(sum_start) = sum_start {
        let add_elapsed = add_start.elapsed();
        if add_elapsed.as_millis() >= 1_000 || done <= 8 || done.is_multiple_of(16) {
            eprintln!(
                "spenso_profile execute.sum_progress done={} total={} add_ms={:.3} elapsed_ms={:.3}",
                done,
                total,
                add_elapsed.as_secs_f64() * 1000.0,
                sum_start.elapsed().as_secs_f64() * 1000.0,
            );
        }
    }
}

#[cfg(feature = "shadowing")]
#[derive(Clone, Copy, Debug, Default)]
struct AtomSumShapeStats {
    logical_entries: usize,
    entries: usize,
    zero_entries: usize,
    top_level_sum_entries: usize,
    total_terms: usize,
    max_terms: usize,
    total_bytes: usize,
    max_bytes: usize,
}

#[cfg(feature = "shadowing")]
impl AtomSumShapeStats {
    fn observe_atom(&mut self, atom: AtomView<'_>, include_bytes: bool) {
        let terms = atom.nterms();
        self.entries += 1;
        self.zero_entries += usize::from(atom.is_zero());
        self.top_level_sum_entries += usize::from(matches!(atom, AtomView::Add(_)));
        self.total_terms += terms;
        self.max_terms = self.max_terms.max(terms);
        if include_bytes {
            let bytes = atom.get_byte_size();
            self.total_bytes += bytes;
            self.max_bytes = self.max_bytes.max(bytes);
        }
    }

    fn merge(&mut self, other: Self) {
        self.logical_entries += other.logical_entries;
        self.entries += other.entries;
        self.zero_entries += other.zero_entries;
        self.top_level_sum_entries += other.top_level_sum_entries;
        self.total_terms += other.total_terms;
        self.max_terms = self.max_terms.max(other.max_terms);
        self.total_bytes += other.total_bytes;
        self.max_bytes = self.max_bytes.max(other.max_bytes);
    }
}

#[cfg(feature = "shadowing")]
trait AtomEntryShape {
    fn observe_atom_shape(self, stats: &mut AtomSumShapeStats, include_bytes: bool);
}

#[cfg(feature = "shadowing")]
impl AtomEntryShape for AtomView<'_> {
    fn observe_atom_shape(self, stats: &mut AtomSumShapeStats, include_bytes: bool) {
        stats.observe_atom(self, include_bytes);
    }
}

#[cfg(feature = "shadowing")]
impl<T> AtomEntryShape for AtomViewOrConcrete<'_, T> {
    fn observe_atom_shape(self, stats: &mut AtomSumShapeStats, include_bytes: bool) {
        if let AtomViewOrConcrete::Atom(atom) = self {
            stats.observe_atom(atom, include_bytes);
        }
    }
}

#[cfg(feature = "shadowing")]
impl<T> AtomEntryShape for &T {
    fn observe_atom_shape(self, _stats: &mut AtomSumShapeStats, _include_bytes: bool) {}
}

#[cfg(feature = "shadowing")]
impl<T> AtomEntryShape for RealOrComplexRef<'_, T> {
    fn observe_atom_shape(self, _stats: &mut AtomSumShapeStats, _include_bytes: bool) {}
}

#[cfg(feature = "shadowing")]
trait AtomSumShapeDiagnostics {
    fn atom_sum_shape_stats(&self, include_bytes: bool) -> AtomSumShapeStats;
}

#[cfg(not(feature = "shadowing"))]
trait AtomSumShapeDiagnostics {}

#[cfg(not(feature = "shadowing"))]
impl<T> AtomSumShapeDiagnostics for T {}

#[cfg(feature = "shadowing")]
impl<T> AtomSumShapeDiagnostics for T
where
    T: HasStructure + IteratableTensor,
    for<'a> T::Data<'a>: AtomEntryShape,
{
    fn atom_sum_shape_stats(&self, include_bytes: bool) -> AtomSumShapeStats {
        let mut stats = AtomSumShapeStats {
            logical_entries: self.structure().size().unwrap_or(0),
            ..Default::default()
        };
        for (_, entry) in self.iter_flat() {
            entry.observe_atom_shape(&mut stats, include_bytes);
        }
        stats
    }
}

#[cfg(feature = "shadowing")]
fn scalar_atom_sum_shape_stats<S: 'static>(
    scalar: &S,
    include_bytes: bool,
) -> Option<AtomSumShapeStats> {
    let atom = (scalar as &dyn std::any::Any).downcast_ref::<Atom>()?;
    let mut stats = AtomSumShapeStats {
        logical_entries: 1,
        ..Default::default()
    };
    stats.observe_atom(atom.as_view(), include_bytes);
    Some(stats)
}

#[cfg(feature = "shadowing")]
fn log_sum_leaf_atom_shapes<K, FK, Aind, Store, LT, L>(
    store: &Store,
    graph: &NetworkGraph<K, FK, Aind>,
    targets: &[(NodeIndex, &NetworkLeaf<K, Aind>)],
    lib: &L,
) where
    Store: NetworkStoreAccess,
    Store::Scalar: 'static,
    Store::Tensor: AtomSumShapeDiagnostics + From<LT::WithIndices> + HasStructure,
    LT: LibraryTensor + Clone,
    L: Library<<Store::Tensor as HasStructure>::Structure, Key = K, Value = PermutedStructure<LT>>,
    K: Display + Debug,
    FK: Debug,
    Aind: AbsInd,
    LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    let include_bytes = profile::atom_bytes();
    let mut leaf_stats = Vec::new();
    for (offset, (node_id, leaf)) in targets.iter().enumerate() {
        let (kind, stats) = match leaf {
            NetworkLeaf::Scalar(index) => (
                "scalar",
                scalar_atom_sum_shape_stats(store.scalar_ref(*index), include_bytes)
                    .unwrap_or_default(),
            ),
            NetworkLeaf::LocalTensor(index) => (
                "tensor",
                store.tensor(*index).atom_sum_shape_stats(include_bytes),
            ),
            NetworkLeaf::TensorSum(indices) => {
                let mut stats = AtomSumShapeStats::default();
                for index in indices {
                    stats.merge(store.tensor(*index).atom_sum_shape_stats(include_bytes));
                }
                ("tensor_sum", stats)
            }
            NetworkLeaf::TensorTerm(term) => {
                let mut stats = store
                    .tensor(term.tensor)
                    .atom_sum_shape_stats(include_bytes);
                if let Some(scalar) = term.scalar {
                    stats.merge(
                        scalar_atom_sum_shape_stats(store.scalar_ref(scalar), include_bytes)
                            .unwrap_or_default(),
                    );
                }
                ("tensor_term", stats)
            }
            NetworkLeaf::TensorTermSum(terms) => {
                let mut stats = AtomSumShapeStats::default();
                for term in terms {
                    stats.merge(
                        store
                            .tensor(term.tensor)
                            .atom_sum_shape_stats(include_bytes),
                    );
                    if let Some(scalar) = term.scalar {
                        stats.merge(
                            scalar_atom_sum_shape_stats(store.scalar_ref(scalar), include_bytes)
                                .unwrap_or_default(),
                        );
                    }
                }
                ("tensor_term_sum", stats)
            }
            NetworkLeaf::LibraryKey { .. } => {
                let tensor = Store::Tensor::from(
                    graph
                        .get_lib_data::<_, LT, L>(lib, *node_id)
                        .expect("library target leaf should resolve to library data"),
                );
                ("library", tensor.atom_sum_shape_stats(include_bytes))
            }
        };
        leaf_stats.push((offset, kind, stats));
    }

    let inspected_atom_leaves = leaf_stats
        .iter()
        .filter(|(_, _, stats)| stats.entries > 0)
        .count();
    let leaves_with_top_level_sums = leaf_stats
        .iter()
        .filter(|(_, _, stats)| stats.top_level_sum_entries > 0)
        .count();
    let atom_entries = leaf_stats
        .iter()
        .map(|(_, _, stats)| stats.entries)
        .sum::<usize>();
    let logical_entries = leaf_stats
        .iter()
        .map(|(_, _, stats)| stats.logical_entries)
        .sum::<usize>();
    let zero_entries = leaf_stats
        .iter()
        .map(|(_, _, stats)| stats.zero_entries)
        .sum::<usize>();
    let structural_zero_entries = logical_entries.saturating_sub(atom_entries);
    let nonzero_entries = atom_entries.saturating_sub(zero_entries);
    let nonzero_density = if logical_entries == 0 {
        0.0
    } else {
        nonzero_entries as f64 / logical_entries as f64
    };
    let top_level_sum_entries = leaf_stats
        .iter()
        .map(|(_, _, stats)| stats.top_level_sum_entries)
        .sum::<usize>();
    let total_terms = leaf_stats
        .iter()
        .map(|(_, _, stats)| stats.total_terms)
        .sum::<usize>();
    let max_terms = leaf_stats
        .iter()
        .map(|(_, _, stats)| stats.max_terms)
        .max()
        .unwrap_or(0);
    let total_bytes = leaf_stats
        .iter()
        .map(|(_, _, stats)| stats.total_bytes)
        .sum::<usize>();
    let max_bytes = leaf_stats
        .iter()
        .map(|(_, _, stats)| stats.max_bytes)
        .max()
        .unwrap_or(0);

    eprintln!(
        "spenso_profile execute.sum_atom_leaves leaves={} inspected_atom_leaves={} leaves_with_top_level_sums={} logical_entries={} atom_entries={} zero_entries={} structural_zero_entries={} nonzero_entries={} nonzero_density={:.6} top_level_sum_entries={} total_terms={} max_terms={} bytes_enabled={} total_bytes={} max_bytes={}",
        targets.len(),
        inspected_atom_leaves,
        leaves_with_top_level_sums,
        logical_entries,
        atom_entries,
        zero_entries,
        structural_zero_entries,
        nonzero_entries,
        nonzero_density,
        top_level_sum_entries,
        total_terms,
        max_terms,
        include_bytes,
        total_bytes,
        max_bytes,
    );

    leaf_stats.sort_by_key(|(_, _, stats)| {
        std::cmp::Reverse((
            stats.total_bytes,
            stats.total_terms,
            stats.entries.saturating_sub(stats.zero_entries),
            stats.top_level_sum_entries,
        ))
    });
    for (rank, (offset, kind, stats)) in leaf_stats
        .into_iter()
        .filter(|(_, _, stats)| stats.entries > 0)
        .take(8)
        .enumerate()
    {
        let structural_zero_entries = stats.logical_entries.saturating_sub(stats.entries);
        let nonzero_entries = stats.entries.saturating_sub(stats.zero_entries);
        let nonzero_density = if stats.logical_entries == 0 {
            0.0
        } else {
            nonzero_entries as f64 / stats.logical_entries as f64
        };
        eprintln!(
            "spenso_profile execute.sum_atom_leaf rank={} leaf={} kind={} logical_entries={} atom_entries={} zero_entries={} structural_zero_entries={} nonzero_entries={} nonzero_density={:.6} top_level_sum_entries={} total_terms={} max_terms={} bytes_enabled={} total_bytes={} max_bytes={}",
            rank,
            offset,
            kind,
            stats.logical_entries,
            stats.entries,
            stats.zero_entries,
            structural_zero_entries,
            nonzero_entries,
            nonzero_density,
            stats.top_level_sum_entries,
            stats.total_terms,
            stats.max_terms,
            include_bytes,
            stats.total_bytes,
            stats.max_bytes,
        );
    }
}

fn balanced_ref_sum<T>(mut terms: Vec<T>, sum_start: Option<std::time::Instant>) -> T
where
    T: Ref + for<'a> AddAssign<<T as Ref>::Ref<'a>>,
{
    debug_assert!(!terms.is_empty());

    let total = terms.len();
    let mut done = 1usize;
    while terms.len() > 1 {
        let mut next = Vec::with_capacity(terms.len().div_ceil(2));
        let mut iter = terms.into_iter();

        while let Some(mut lhs) = iter.next() {
            if let Some(rhs) = iter.next() {
                let add_start = std::time::Instant::now();
                lhs += rhs.refer();
                done += 1;
                log_sum_add(sum_start, done, total, add_start);
            }
            next.push(lhs);
        }

        terms = next;
    }

    terms.pop().expect("balanced sum has at least one term")
}

fn try_balanced_scalar_sum<K, Aind, Store>(
    store: &mut Store,
    targets: &[(NodeIndex, &NetworkLeaf<K, Aind>)],
    sum_start: Option<std::time::Instant>,
) -> Option<NetworkLeaf<K, Aind>>
where
    Store: NetworkStoreAccess,
    Store::Scalar: Clone + Ref + for<'a> AddAssign<<Store::Scalar as Ref>::Ref<'a>>,
{
    let mut terms = Vec::with_capacity(targets.len());
    for (_, leaf) in targets {
        let NetworkLeaf::Scalar(index) = leaf else {
            return None;
        };
        terms.push(store.scalar_ref(*index).clone());
    }

    let result = balanced_ref_sum(terms, sum_start);
    Some(NetworkLeaf::Scalar(store.push_scalar(result).into()))
}

fn try_balanced_tensor_sum<K, FK, Aind, Store, LT, L>(
    store: &mut Store,
    graph: &NetworkGraph<K, FK, Aind>,
    targets: &[(NodeIndex, &NetworkLeaf<K, Aind>)],
    lib: &L,
    sum_start: Option<std::time::Instant>,
) -> Option<NetworkLeaf<K, Aind>>
where
    Store: NetworkStoreAccess,
    Store::Tensor: HasStructure
        + TensorStructure
        + Clone
        + Ref
        + FastTensorSum
        + From<LT::WithIndices>
        + ScalarMul<Store::Scalar, Output = Store::Tensor>
        + for<'a> AddAssign<<Store::Tensor as Ref>::Ref<'a>>,
    Store::Scalar: Clone,
    LT: LibraryTensor + Clone,
    L: Library<<Store::Tensor as HasStructure>::Structure, Key = K, Value = PermutedStructure<LT>>,
    K: Display + Debug,
    FK: Debug,
    Aind: AbsInd,
    LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    let mut terms = Vec::with_capacity(targets.len());
    for (node_id, leaf) in targets {
        match leaf {
            NetworkLeaf::Scalar(_) => return None,
            NetworkLeaf::LocalTensor(index) => {
                terms.push(TensorTerm::tensor(*index));
            }
            NetworkLeaf::TensorSum(indices) => {
                for index in indices {
                    terms.push(TensorTerm::tensor(*index));
                }
            }
            NetworkLeaf::TensorTerm(term) => terms.push(term.clone()),
            NetworkLeaf::TensorTermSum(indices) => terms.extend(indices.iter().cloned()),
            NetworkLeaf::LibraryKey { .. } => {
                let tensor = Store::Tensor::from(graph.get_lib_data::<_, LT, L>(lib, *node_id)?);
                terms.push(TensorTerm::tensor(store.push_tensor(tensor)));
            }
        }
    }

    let all_scalar_terms = terms
        .iter()
        .all(|term| store.tensor(term.tensor).scalar_ref().is_some());
    let keep_lazy = profile::lazy_tensor_sums()
        && !all_scalar_terms
        && terms.len() >= MIN_LAZY_TENSOR_SUM_TERMS;

    if profile::enabled()
        && (keep_lazy || targets.len() > LARGE_SUM_PROFILE_THRESHOLD || profile::verbose())
    {
        eprintln!(
            "spenso_profile execute.sum_tensor_mode leaves={} terms={} mode={}",
            targets.len(),
            terms.len(),
            if keep_lazy { "lazy" } else { "materialized" },
        );
    }

    if !keep_lazy {
        Some(materialize_tensor_terms(store, &terms, sum_start))
    } else {
        Some(tensor_term_sum_leaf(&terms))
    }
}

fn tensor_term_sum_leaf<K, Aind>(terms: &[TensorTerm]) -> NetworkLeaf<K, Aind> {
    debug_assert!(!terms.is_empty());

    if terms.len() == 1 {
        let term = terms[0].clone();
        return match term.scalar {
            Some(_) => NetworkLeaf::TensorTerm(term),
            None => NetworkLeaf::LocalTensor(term.tensor),
        };
    }

    if terms.iter().all(|term| term.scalar.is_none()) {
        NetworkLeaf::TensorSum(terms.iter().map(|term| term.tensor).collect())
    } else {
        NetworkLeaf::TensorTermSum(terms.to_vec())
    }
}

fn materialize_tensor_terms<K, Aind, Store>(
    store: &mut Store,
    terms: &[TensorTerm],
    sum_start: Option<std::time::Instant>,
) -> NetworkLeaf<K, Aind>
where
    Store: NetworkStoreAccess,
    Store::Tensor: Clone
        + Ref
        + FastTensorSum
        + ScalarMul<Store::Scalar, Output = Store::Tensor>
        + for<'a> AddAssign<<Store::Tensor as Ref>::Ref<'a>>,
    Store::Scalar: Clone,
{
    if terms.iter().all(|term| term.scalar.is_none()) {
        let indices = terms.iter().map(|term| term.tensor).collect::<Vec<_>>();
        return materialize_tensor_sum(store, &indices, sum_start);
    }

    let mut materialized_terms = Vec::with_capacity(terms.len());
    for term in terms {
        let tensor = match term.scalar {
            Some(scalar) => store
                .tensor(term.tensor)
                .scalar_mul(store.scalar_ref(scalar))
                .expect("scaled tensor term should support scalar multiplication"),
            None => store.tensor(term.tensor).clone(),
        };
        materialized_terms.push(store.push_tensor(tensor));
    }

    materialize_tensor_sum(store, &materialized_terms, sum_start)
}

fn materialize_tensor_sum<K, Aind, Store>(
    store: &mut Store,
    indices: &[usize],
    sum_start: Option<std::time::Instant>,
) -> NetworkLeaf<K, Aind>
where
    Store: NetworkStoreAccess,
    Store::Tensor: Clone + Ref + FastTensorSum + for<'a> AddAssign<<Store::Tensor as Ref>::Ref<'a>>,
{
    if let Some(result) = {
        let terms = indices
            .iter()
            .map(|index| store.tensor(*index))
            .collect::<Vec<_>>();
        Store::Tensor::fast_tensor_sum(&terms, sum_start)
    } {
        return NetworkLeaf::LocalTensor(store.push_tensor(result));
    }

    let terms = indices
        .iter()
        .map(|index| store.tensor(*index).clone())
        .collect::<Vec<_>>();
    let result = balanced_ref_sum(terms, sum_start);
    NetworkLeaf::LocalTensor(store.push_tensor(result))
}

#[cfg(feature = "shadowing")]
fn try_atom_scalar_sum<K, Aind, Store>(
    store: &mut Store,
    targets: &[(NodeIndex, &NetworkLeaf<K, Aind>)],
) -> Option<NetworkLeaf<K, Aind>>
where
    Store: NetworkStoreAccess,
    Store::Scalar: 'static,
{
    use std::{
        any::{Any, TypeId},
        fs::File,
        io::BufWriter,
    };

    use symbolica::{
        atom::Atom,
        streaming::{TermStreamer, TermStreamerConfig},
    };

    if TypeId::of::<Store::Scalar>() != TypeId::of::<Atom>() {
        return None;
    }

    const MIN_ATOM_ADD_MANY_LEAVES: usize = 4;
    const MIN_STREAMED_SUM_LEAVES: usize = 32;
    if targets.len() < MIN_ATOM_ADD_MANY_LEAVES {
        return None;
    }

    let atoms = targets
        .iter()
        .map(|(_, leaf)| {
            let NetworkLeaf::Scalar(index) = leaf else {
                return None;
            };
            (store.scalar_ref(*index) as &dyn Any)
                .downcast_ref::<Atom>()
                .cloned()
        })
        .collect::<Option<Vec<_>>>()?;

    let start = profile::enabled().then(std::time::Instant::now);

    let result = if targets.len() >= MIN_STREAMED_SUM_LEAVES {
        if profile::enabled() {
            eprintln!(
                "spenso_profile execute.sum_term_stream_start leaves={}",
                targets.len(),
            );
        }

        let mut stream = TermStreamer::<BufWriter<File>>::new(
            TermStreamerConfig::new()
                .path(std::env::temp_dir().to_string_lossy().into_owned())
                .max_mem_bytes(4 * 1024 * 1024 * 1024),
        );

        for atom in atoms {
            stream.push(atom);
        }

        let result = stream.to_expression();
        if let Some(start) = start {
            eprintln!(
                "spenso_profile execute.sum_term_stream_done leaves={} terms={} bytes={} elapsed_ms={:.3}",
                targets.len(),
                result.nterms(),
                result.as_view().get_byte_size(),
                start.elapsed().as_secs_f64() * 1000.0,
            );
        }
        result
    } else {
        if profile::verbose() {
            eprintln!(
                "spenso_profile execute.sum_atom_add_many_start leaves={}",
                targets.len(),
            );
        }

        let result = Atom::add_many(&atoms);
        if profile::verbose()
            && let Some(start) = start
        {
            eprintln!(
                "spenso_profile execute.sum_atom_add_many_done leaves={} terms={} bytes={} elapsed_ms={:.3}",
                targets.len(),
                result.nterms(),
                result.as_view().get_byte_size(),
                start.elapsed().as_secs_f64() * 1000.0,
            );
        }
        result
    };

    let boxed: Box<dyn Any> = Box::new(result);
    let scalar = *boxed.downcast::<Store::Scalar>().ok()?;
    Some(NetworkLeaf::Scalar(store.push_scalar(scalar).into()))
}

pub struct Sequential;
pub struct Parallel;
pub struct SequentialRef;
pub struct SequentialExtract;

pub enum ExecutionMode {
    Sequential,
    Parallel,
    SequentialRef,
    SequentialExtract,
}

pub struct Steps<const N: usize> {}
pub struct StepsDebug<const N: usize> {}

impl<const N: usize, E, L, FL, FK: Debug, K: Debug, Aind: AbsInd>
    ExecutionStrategy<E, FL, L, K, FK, Aind> for StepsDebug<N>
where
    E: ExecuteOp<FL, L, K, FK, Aind>,
    K: Clone,
    FK: Clone,
{
    fn execute_all<C: ContractionStrategy<E, L, K, FK, Aind>>(
        executor: &mut E,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        fnlib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        if <C as ContractionStrategy<E, L, K, FK, Aind>>::SUPPORTS_PARTIAL_GRAPH_REWRITE {
            let mut ignored: SuBitGraph = graph.graph.empty_subgraph();
            for _ in 0..N {
                while executor.execute_self_loop_traces_ignoring(graph, lib, &mut ignored)? {}

                profile::bump(Counter::ExecuteIteration, 1);
                let planned = plan_ready_operation_batch(graph, &ignored);
                if planned.is_empty() {
                    break;
                }

                let mut did_progress = false;
                for operation in planned {
                    if execute_operation_in_place::<C, E, FL, L, K, FK, Aind>(
                        executor,
                        graph,
                        &operation,
                        lib,
                        fnlib,
                        &mut ignored,
                    )? {
                        did_progress = true;
                        break;
                    }
                }

                if !did_progress {
                    break;
                }
            }

            if !ignored.is_empty() {
                graph.delete(&ignored);
            }

            return Ok(());
        }

        for _ in 0..N {
            while executor.execute_self_loop_traces(graph, lib)? {}

            // find the *one* ready op
            if let Some((mut extracted_graph, _op)) = graph.extract_next_ready_op() {
                println!(
                    "Extracted_graph: {}",
                    extracted_graph.dot_impl(
                        |s| s.to_string(),
                        |_| "".to_string(),
                        |s| s.to_string(),
                        |_| "".to_string()
                    )
                );
                println!(
                    "Graph: {}",
                    graph.dot_impl(
                        |s| s.to_string(),
                        |_| "".to_string(),
                        |s| s.to_string(),
                        |_| "".to_string()
                    )
                );
                // execute + splice
                let operation =
                    first_executable_operation(&mut extracted_graph).ok_or_else(|| {
                        TensorNetworkError::Other(eyre!(
                            "extracted graph has no executable operation"
                        ))
                    })?;
                let replacement =
                    executor.execute::<C>(&extracted_graph, &operation, lib, fnlib)?;
                collapse_operation_subgraph(&mut extracted_graph, &operation, replacement)?;
                println!(
                    "Replacement Graph: {}",
                    extracted_graph.dot_impl(
                        |s| s.to_string(),
                        |_| "".to_string(),
                        |s| s.to_string(),
                        |_| "".to_string()
                    )
                );

                graph.splice_descendents_of(extracted_graph);
            }
        }

        Ok(())
    }
}

impl<const N: usize, E, FL, L, K, FK, Aind: AbsInd> ExecutionStrategy<E, FL, L, K, FK, Aind>
    for Steps<N>
where
    E: ExecuteOp<FL, L, K, FK, Aind>,
    K: Clone + Debug,
    FK: Clone + Debug,
{
    fn execute_all<C: ContractionStrategy<E, L, K, FK, Aind>>(
        executor: &mut E,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        fnlib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        if <C as ContractionStrategy<E, L, K, FK, Aind>>::SUPPORTS_PARTIAL_GRAPH_REWRITE {
            let mut ignored: SuBitGraph = graph.graph.empty_subgraph();
            for _ in 0..N {
                while executor.execute_self_loop_traces_ignoring(graph, lib, &mut ignored)? {}

                profile::bump(Counter::ExecuteIteration, 1);
                let planned = plan_ready_operation_batch(graph, &ignored);
                if planned.is_empty() {
                    break;
                }

                let mut did_progress = false;
                for operation in planned {
                    if execute_operation_in_place::<C, E, FL, L, K, FK, Aind>(
                        executor,
                        graph,
                        &operation,
                        lib,
                        fnlib,
                        &mut ignored,
                    )? {
                        did_progress = true;
                        break;
                    }
                }

                if !did_progress {
                    break;
                }
            }

            if !ignored.is_empty() {
                graph.delete(&ignored);
            }

            return Ok(());
        }

        for _ in 0..N {
            while executor.execute_self_loop_traces(graph, lib)? {}

            // find the *one* ready op
            if let Some((mut extracted_graph, _op)) = graph.extract_next_ready_op() {
                // execute + splice
                let operation =
                    first_executable_operation(&mut extracted_graph).ok_or_else(|| {
                        TensorNetworkError::Other(eyre!(
                            "extracted graph has no executable operation"
                        ))
                    })?;
                let replacement =
                    executor.execute::<C>(&extracted_graph, &operation, lib, fnlib)?;
                collapse_operation_subgraph(&mut extracted_graph, &operation, replacement)?;
                graph.splice_descendents_of(extracted_graph);
            }
        }

        Ok(())
    }
}

impl<E, L, FL, K, FK, Aind: AbsInd> ExecutionStrategy<E, FL, L, K, FK, Aind> for SequentialExtract
where
    E: ExecuteOp<FL, L, K, FK, Aind>,
    FK: Clone + Debug,
    K: Clone + Debug,
{
    fn execute_all<C: ContractionStrategy<E, L, K, FK, Aind>>(
        executor: &mut E,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        fnlib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        while {
            // find the *one* ready op
            profile::bump(Counter::ExecuteIteration, 1);
            let ready = {
                let _span = profile::span(Timer::ExecuteFindReady);
                graph.extract_next_ready_op()
            };
            if let Some((mut extracted_graph, _op)) = ready {
                profile::bump(Counter::ExecuteFound, 1);
                // execute + splice
                let operation =
                    first_executable_operation(&mut extracted_graph).ok_or_else(|| {
                        TensorNetworkError::Other(eyre!(
                            "extracted graph has no executable operation"
                        ))
                    })?;
                let replacement =
                    executor.execute::<C>(&extracted_graph, &operation, lib, fnlib)?;
                collapse_operation_subgraph(&mut extracted_graph, &operation, replacement)?;
                graph.splice_descendents_of(extracted_graph);
                true
            } else {
                false
            }
        } {}

        Ok(())
    }
}

impl<E, L, FL, K, FK, Aind: AbsInd> ExecutionStrategy<E, FL, L, K, FK, Aind> for SequentialRef
where
    E: ExecuteOp<FL, L, K, FK, Aind>,
    FK: Clone + Debug,
    K: Clone + Debug,
{
    fn execute_all<C: ContractionStrategy<E, L, K, FK, Aind>>(
        executor: &mut E,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        fnlib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        while {
            profile::bump(Counter::ExecuteIteration, 1);
            let ready = {
                let _span = profile::span(Timer::ExecuteFindReady);
                graph.extract_next_ready_ref_op()
            };
            if let Some((mut extracted_graph, _op)) = ready {
                profile::bump(Counter::ExecuteFound, 1);
                let operation =
                    first_executable_operation(&mut extracted_graph).ok_or_else(|| {
                        TensorNetworkError::Other(eyre!(
                            "extracted graph has no executable operation"
                        ))
                    })?;
                let replacement =
                    executor.execute::<C>(&extracted_graph, &operation, lib, fnlib)?;
                collapse_operation_subgraph(&mut extracted_graph, &operation, replacement)?;
                graph.splice_descendents_of(extracted_graph);
                true
            } else {
                false
            }
        } {}

        Ok(())
    }
}

impl<E, L, FL, K, FK, Aind: AbsInd> ExecutionStrategy<E, FL, L, K, FK, Aind> for Sequential
where
    E: ExecuteOp<FL, L, K, FK, Aind>,
    FK: Clone + Debug,
    K: Clone + Debug,
{
    fn execute_all<C: ContractionStrategy<E, L, K, FK, Aind>>(
        executor: &mut E,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        fnlib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        let mut ignored: SuBitGraph = graph.graph.empty_subgraph();

        if <C as ContractionStrategy<E, L, K, FK, Aind>>::SUPPORTS_PARTIAL_GRAPH_REWRITE {
            loop {
                if executor.execute_self_loop_traces_ignoring(graph, lib, &mut ignored)? {
                    continue;
                }

                profile::bump(Counter::ExecuteIteration, 1);
                let planned = plan_ready_operation_batch(graph, &ignored);
                if planned.is_empty() {
                    break;
                }

                let mut did_progress = false;
                for operation in planned {
                    if execute_operation_in_place::<C, E, FL, L, K, FK, Aind>(
                        executor,
                        graph,
                        &operation,
                        lib,
                        fnlib,
                        &mut ignored,
                    )? {
                        did_progress = true;
                        break;
                    }
                }

                if !did_progress {
                    break;
                }
            }

            if !ignored.is_empty() {
                graph.delete(&ignored);
            }

            return Ok(());
        }

        let mut batch_index = 0usize;

        loop {
            if executor.execute_self_loop_traces_ignoring(graph, lib, &mut ignored)? {
                continue;
            }

            profile::bump(Counter::ExecuteIteration, 1);
            let planned = plan_ready_operation_batch(graph, &ignored);

            if planned.is_empty() {
                break;
            }

            profile::bump(Counter::ExecuteFound, planned.len() as u64);
            let mut replacements = Vec::with_capacity(planned.len());
            let profile_batch = profile::enabled();
            let batch_start = if profile_batch {
                Some(std::time::Instant::now())
            } else {
                None
            };
            for (op_index, planned_op) in planned.iter().enumerate() {
                let op_start = if profile_batch {
                    if planned.len() <= 1024 && (op_index < 16 || op_index % 64 == 0) {
                        eprintln!(
                            "spenso_profile execute.batch_op_start batch={} op_index={} total={} op={} leaves={} subgraph_hedges={}",
                            batch_index,
                            op_index,
                            planned.len(),
                            planned_op.op().display_with(|fun| fun.to_string()),
                            planned_op.leaf_count(),
                            planned_op.subgraph().n_included(),
                        );
                    }
                    Some(std::time::Instant::now())
                } else {
                    None
                };
                replacements.push(executor.execute::<C>(graph, planned_op, lib, fnlib)?);
                if let Some(op_start) = op_start {
                    let elapsed = op_start.elapsed();
                    if profile::verbose() || elapsed.as_millis() >= 100 {
                        eprintln!(
                            "spenso_profile execute.slow_op batch={} op_index={} elapsed_ms={:.3} op={} leaves={} subgraph_hedges={}",
                            batch_index,
                            op_index,
                            elapsed.as_secs_f64() * 1000.0,
                            planned_op.op().display_with(|fun| fun.to_string()),
                            planned_op.leaf_count(),
                            planned_op.subgraph().n_included(),
                        );
                    }
                    if planned.len() <= 1024 && (op_index + 1) % 64 == 0 {
                        eprintln!(
                            "spenso_profile execute.batch_progress batch={} done={} total={} elapsed_ms={:.3}",
                            batch_index,
                            op_index + 1,
                            planned.len(),
                            batch_start
                                .map(|start| start.elapsed().as_secs_f64() * 1000.0)
                                .unwrap_or_default(),
                        );
                    }
                }
            }
            if let Some(batch_start) = batch_start {
                eprintln!(
                    "spenso_profile execute.batch_done batch={} ops={} elapsed_ms={:.3}",
                    batch_index,
                    planned.len(),
                    batch_start.elapsed().as_secs_f64() * 1000.0,
                );
            }

            for (planned_op, replacement) in planned.into_iter().zip(replacements) {
                graph
                    .identify_subgraph_nodes_without_deleting_self_edges(
                        planned_op.subgraph(),
                        NetworkNode::Leaf(replacement),
                        &mut ignored,
                    )
                    .ok_or_else(|| {
                        TensorNetworkError::Other(eyre!(
                            "ready operation subgraph did not contain any nodes"
                        ))
                    })?;
            }
            graph.finish_deferred_node_identifications();
            batch_index += 1;
        }

        if !ignored.is_empty() {
            graph.delete(&ignored);
        }

        Ok(())
    }
}

struct ParallelExecutionOutput<T, Sc, K, Aind> {
    replacement: NetworkLeaf<K, Aind>,
    tensors: Vec<T>,
    scalars: Vec<Sc>,
    scalar_aliases: Vec<Option<Sc>>,
}

fn remap_parallel_replacement<K, Aind>(
    replacement: &mut NetworkLeaf<K, Aind>,
    base_tensors: usize,
    base_scalars: usize,
    tensor_offset: usize,
    scalar_offset: usize,
) {
    fn rebase_index(index: &mut usize, base: usize, offset: usize) {
        if *index >= base {
            *index = offset + (*index - base);
        }
    }

    replacement.map_scalar_refs(|scalar| scalar.rebase_from(base_scalars, scalar_offset));

    match replacement {
        NetworkLeaf::LocalTensor(index) => rebase_index(index, base_tensors, tensor_offset),
        NetworkLeaf::TensorSum(indices) => {
            for index in indices {
                rebase_index(index, base_tensors, tensor_offset);
            }
        }
        NetworkLeaf::TensorTerm(term) => {
            rebase_index(&mut term.tensor, base_tensors, tensor_offset);
        }
        NetworkLeaf::TensorTermSum(terms) => {
            for term in terms {
                rebase_index(&mut term.tensor, base_tensors, tensor_offset);
            }
        }
        NetworkLeaf::Scalar(_) | NetworkLeaf::LibraryKey { .. } => {}
    }
}

impl Parallel {
    #[allow(clippy::result_large_err)]
    pub fn execute_all<T, Sc, L, FL, K, FK, Aind, C>(
        executor: &mut NetworkStore<T, Sc>,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        fnlib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        NetworkStore<T, Sc>: ExecuteOp<FL, L, K, FK, Aind>,
        for<'a> NetworkStoreOverlay<'a, T, Sc>: ExecuteOp<FL, L, K, FK, Aind>,
        C: ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind>,
        for<'a> C: ContractionStrategy<NetworkStoreOverlay<'a, T, Sc>, L, K, FK, Aind>,
        T: Clone + Send + Sync,
        Sc: Clone + Send + Sync,
        L: Sync,
        FL: Sync,
        K: Clone + Debug + Display + Send + Sync,
        FK: Clone + Debug + Display + Send + Sync,
        Aind: AbsInd + Send + Sync,
    {
        if <C as ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind>>::SUPPORTS_PARTIAL_GRAPH_REWRITE
        {
            return <Sequential as ExecutionStrategy<
                NetworkStore<T, Sc>,
                FL,
                L,
                K,
                FK,
                Aind,
            >>::execute_all::<C>(executor, graph, lib, fnlib);
        }

        let mut ignored: SuBitGraph = graph.graph.empty_subgraph();

        loop {
            profile::bump(Counter::ExecuteIteration, 1);
            let planned = plan_ready_operation_batch(graph, &ignored);

            if planned.is_empty() {
                break;
            }

            profile::bump(Counter::ExecuteFound, planned.len() as u64);
            let replacements = if planned.len() == 1 {
                vec![executor.execute::<C>(graph, &planned[0], lib, fnlib)?]
            } else {
                let base_tensors = executor.tensors.len();
                let base_scalars = executor.scalar.len();
                let base_executor = &*executor;
                let worker_ops = planned.clone();

                let outputs = worker_ops
                    .into_par_iter()
                    .map(|planned_op| -> Result<_, Box<TensorNetworkError<K, FK>>> {
                        let mut local = NetworkStoreOverlay::new(base_executor);
                        let replacement = local
                            .execute::<C>(graph, &planned_op, lib, fnlib)
                            .map_err(Box::new)?;
                        let (tensors, scalars, scalar_aliases) = local.into_additions();

                        Ok(ParallelExecutionOutput {
                            replacement,
                            tensors,
                            scalars,
                            scalar_aliases,
                        })
                    })
                    .collect::<Result<Vec<_>, Box<TensorNetworkError<K, FK>>>>()
                    .map_err(|err| *err)?;

                let mut replacements = Vec::with_capacity(outputs.len());
                for output in outputs {
                    let tensor_offset = executor.tensors.len();
                    let scalar_offset = executor.scalar.len();
                    let mut replacement = output.replacement;
                    remap_parallel_replacement(
                        &mut replacement,
                        base_tensors,
                        base_scalars,
                        tensor_offset,
                        scalar_offset,
                    );

                    executor.tensors.extend(output.tensors);
                    executor.scalar.extend(output.scalars);
                    executor.scalar_aliases.extend(output.scalar_aliases);
                    replacements.push(replacement);
                }

                replacements
            };

            for (planned_op, replacement) in planned.into_iter().zip(replacements) {
                graph
                    .identify_subgraph_nodes_without_deleting_self_edges(
                        planned_op.subgraph(),
                        NetworkNode::Leaf(replacement),
                        &mut ignored,
                    )
                    .ok_or_else(|| {
                        TensorNetworkError::Other(eyre!(
                            "ready operation subgraph did not contain any nodes"
                        ))
                    })?;
            }
            graph.finish_deferred_node_identifications();
        }

        if !ignored.is_empty() {
            graph.delete(&ignored);
        }

        Ok(())
    }
}

pub trait ExecuteOp<FL, L, K, FK, Aind>: Sized {
    // type LibStruct;
    #[allow(clippy::result_large_err)]
    fn execute_self_loop_traces(
        &mut self,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display;

    #[allow(clippy::result_large_err)]
    fn execute_self_loop_traces_ignoring(
        &mut self,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        ignored: &mut SuBitGraph,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display;

    #[allow(clippy::result_large_err)]
    fn execute<C: ContractionStrategy<Self, L, K, FK, Aind>>(
        &mut self,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
        fn_lib: &FL,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display;
}

impl<S, Store: TensorScalarStore, K, FK, Aind: AbsInd> Network<Store, K, FK, Aind>
where
    Store::Tensor: HasStructure<Structure = S>,
{
    #[allow(clippy::result_large_err)]
    pub fn execute<
        Strat: ExecutionStrategy<Store, FL, L, K, FK, Aind>,
        C: ContractionStrategy<Store, L, K, FK, Aind>,
        LT,
        L,
        FL,
    >(
        &mut self,
        lib: &L,
        fn_lib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display + Clone + Debug,
        FK: Display + Clone + Debug,
        L: Library<S, Key = K, Value = PermutedStructure<LT>> + Sync,
        FL: FunctionLibrary<Store::Tensor, Store::Scalar, Key = FK>,
        LT: LibraryTensor<WithIndices = Store::Tensor>,
        Store: ExecuteOp<FL, L, K, FK, Aind>,
    {
        {
            let _span = profile::span(Timer::MergeOps);
            profile::bump(Counter::MergeOps, 1);
            self.merge_ops();
        }
        self.store.execute_self_loop_traces(&mut self.graph, lib)?;
        Strat::execute_all::<C>(&mut self.store, &mut self.graph, lib, fn_lib)?;
        self.state = self.graph.state();
        Ok(())
    }
}

fn execute_next_self_loop_trace<LT, L, K, FK, Aind, Store>(
    store: &mut Store,
    graph: &mut NetworkGraph<K, FK, Aind>,
    lib: &L,
) -> Result<bool, TensorNetworkError<K, FK>>
where
    Store: NetworkStoreAccess,
    Store::Tensor: HasStructure
        + TensorStructure
        + Clone
        + Trace
        + Ref
        + FastTensorSum
        + ScalarMul<Store::Scalar, Output = Store::Tensor>
        + for<'a> AddAssign<<Store::Tensor as Ref>::Ref<'a>>
        + From<LT::WithIndices>,
    Store::Scalar: Clone,
    LT: LibraryTensor + Clone,
    L: Library<<Store::Tensor as HasStructure>::Structure, Key = K, Value = PermutedStructure<LT>>,
    K: Display + Debug,
    FK: Display + Debug,
    Aind: AbsInd,
    LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    let Some(node) = graph.graph.iter_nodes().find_map(|(node, crown, data)| {
        if !matches!(data, NetworkNode::Leaf(_)) {
            return None;
        }

        crown
            .clone()
            .any(|hedge| graph.graph[[&hedge]].is_slot() && graph.graph.is_self_loop(hedge))
            .then_some(node)
    }) else {
        return Ok(false);
    };

    let traced = match &graph.graph[node] {
        NetworkNode::Leaf(NetworkLeaf::LocalTensor(local)) => {
            store.tensor(*local).internal_contract()
        }
        NetworkNode::Leaf(NetworkLeaf::LibraryKey { .. }) => Store::Tensor::from(
            graph
                .get_lib_data(lib, node)
                .expect("library node selected for trace must have library data"),
        )
        .internal_contract(),
        NetworkNode::Leaf(NetworkLeaf::TensorSum(indices)) => {
            let materialized = materialize_tensor_sum::<K, Aind, Store>(store, indices, None);
            let NetworkLeaf::LocalTensor(local) = materialized else {
                unreachable!("materialized tensor sum is a local tensor")
            };
            store.tensor(local).internal_contract()
        }
        NetworkNode::Leaf(NetworkLeaf::TensorTerm(term)) => {
            let materialized =
                materialize_tensor_terms::<K, Aind, Store>(store, std::slice::from_ref(term), None);
            let NetworkLeaf::LocalTensor(local) = materialized else {
                unreachable!("materialized tensor term is a local tensor")
            };
            store.tensor(local).internal_contract()
        }
        NetworkNode::Leaf(NetworkLeaf::TensorTermSum(terms)) => {
            let materialized = materialize_tensor_terms::<K, Aind, Store>(store, terms, None);
            let NetworkLeaf::LocalTensor(local) = materialized else {
                unreachable!("materialized tensor terms are a local tensor")
            };
            store.tensor(local).internal_contract()
        }
        NetworkNode::Leaf(NetworkLeaf::Scalar(_)) => return Ok(false),
        NetworkNode::Op(_) => unreachable!("self-loop trace search only selects tensor leaves"),
    };

    let traced_position = store.push_tensor(traced);
    graph.replace_node_deleting_self_loop_slots(
        node,
        NetworkNode::Leaf(NetworkLeaf::LocalTensor(traced_position)),
    );

    Ok(true)
}

fn execute_next_self_loop_trace_ignoring<LT, L, K, FK, Aind, Store>(
    store: &mut Store,
    graph: &mut NetworkGraph<K, FK, Aind>,
    lib: &L,
    ignored: &mut SuBitGraph,
) -> Result<bool, TensorNetworkError<K, FK>>
where
    Store: NetworkStoreAccess,
    Store::Tensor: HasStructure
        + TensorStructure
        + Clone
        + Trace
        + Ref
        + FastTensorSum
        + ScalarMul<Store::Scalar, Output = Store::Tensor>
        + for<'a> AddAssign<<Store::Tensor as Ref>::Ref<'a>>
        + From<LT::WithIndices>,
    Store::Scalar: Clone,
    LT: LibraryTensor + Clone,
    L: Library<<Store::Tensor as HasStructure>::Structure, Key = K, Value = PermutedStructure<LT>>,
    K: Display + Debug,
    FK: Display + Debug,
    Aind: AbsInd,
    LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    let Some(node) = graph.graph.iter_nodes().find_map(|(node, crown, data)| {
        if !matches!(data, NetworkNode::Leaf(_)) {
            return None;
        }

        crown
            .clone()
            .any(|hedge| {
                !ignored.includes(&hedge)
                    && graph.graph[[&hedge]].is_slot()
                    && graph.graph.is_self_loop(hedge)
            })
            .then_some(node)
    }) else {
        return Ok(false);
    };

    let traced = match &graph.graph[node] {
        NetworkNode::Leaf(NetworkLeaf::LocalTensor(local)) => {
            store.tensor(*local).internal_contract()
        }
        NetworkNode::Leaf(NetworkLeaf::LibraryKey { .. }) => Store::Tensor::from(
            graph
                .get_lib_data(lib, node)
                .expect("library node selected for trace must have library data"),
        )
        .internal_contract(),
        NetworkNode::Leaf(NetworkLeaf::TensorSum(indices)) => {
            let materialized = materialize_tensor_sum::<K, Aind, Store>(store, indices, None);
            let NetworkLeaf::LocalTensor(local) = materialized else {
                unreachable!("materialized tensor sum is a local tensor")
            };
            store.tensor(local).internal_contract()
        }
        NetworkNode::Leaf(NetworkLeaf::TensorTerm(term)) => {
            let materialized =
                materialize_tensor_terms::<K, Aind, Store>(store, std::slice::from_ref(term), None);
            let NetworkLeaf::LocalTensor(local) = materialized else {
                unreachable!("materialized tensor term is a local tensor")
            };
            store.tensor(local).internal_contract()
        }
        NetworkNode::Leaf(NetworkLeaf::TensorTermSum(terms)) => {
            let materialized = materialize_tensor_terms::<K, Aind, Store>(store, terms, None);
            let NetworkLeaf::LocalTensor(local) = materialized else {
                unreachable!("materialized tensor terms are a local tensor")
            };
            store.tensor(local).internal_contract()
        }
        NetworkNode::Leaf(NetworkLeaf::Scalar(_)) => return Ok(false),
        NetworkNode::Op(_) => unreachable!("self-loop trace search only selects tensor leaves"),
    };

    let traced_position = store.push_tensor(traced);
    graph.replace_node_ignoring_self_loop_slots(
        node,
        NetworkNode::Leaf(NetworkLeaf::LocalTensor(traced_position)),
        ignored,
    );

    Ok(true)
}

impl<S, T, Sc, K, FK, Aind: AbsInd> Network<NetworkStore<T, Sc>, K, FK, Aind>
where
    T: HasStructure<Structure = S>,
{
    #[allow(clippy::result_large_err)]
    pub fn execute_parallel<C, LT, L, FL>(
        &mut self,
        lib: &L,
        fn_lib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display + Clone + Debug + Send + Sync,
        FK: Display + Clone + Debug + Send + Sync,
        L: Library<S, Key = K, Value = PermutedStructure<LT>> + Sync,
        FL: FunctionLibrary<T, Sc, Key = FK> + Sync,
        LT: LibraryTensor<WithIndices = T>,
        NetworkStore<T, Sc>: ExecuteOp<FL, L, K, FK, Aind>,
        for<'a> NetworkStoreOverlay<'a, T, Sc>: ExecuteOp<FL, L, K, FK, Aind>,
        C: ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind>,
        for<'a> C: ContractionStrategy<NetworkStoreOverlay<'a, T, Sc>, L, K, FK, Aind>,
        T: Clone + Send + Sync,
        Sc: Clone + Send + Sync,
        Aind: Send + Sync,
    {
        {
            let _span = profile::span(Timer::MergeOps);
            profile::bump(Counter::MergeOps, 1);
            self.merge_ops();
        }
        Parallel::execute_all::<T, Sc, L, FL, K, FK, Aind, C>(
            &mut self.store,
            &mut self.graph,
            lib,
            fn_lib,
        )
    }
}

impl<LT, L, K, FK, FL, Aind, Store> ExecuteOp<FL, L, K, FK, Aind> for Store
where
    LT: LibraryTensor + Clone,
    Store: NetworkStoreAccess,
    Store::Tensor: HasStructure
        + TensorStructure
        + Neg<Output = Store::Tensor>
        + Clone
        + Trace
        + Ref
        + FastTensorSum
        + Contract<LCM = Store::Tensor>
        + ScalarMul<Store::Scalar, Output = Store::Tensor>
        + for<'a> AddAssign<<Store::Tensor as Ref>::Ref<'a>>
        + for<'a> AddAssign<LT::WithIndices>
        + From<LT::WithIndices>
        + AtomSumShapeDiagnostics,
    L: Library<<Store::Tensor as HasStructure>::Structure, Key = K, Value = PermutedStructure<LT>>,
    Store::Scalar: Neg<Output = Store::Scalar>
        + RefOne
        + Div<Output = Store::Scalar>
        + for<'a> AddAssign<<Store::Scalar as Ref>::Ref<'a>>
        + Clone
        + for<'a> AddAssign<<Store::Tensor as HasStructure>::ScalarRef<'a>>
        + From<<Store::Tensor as HasStructure>::Scalar>
        + Ref
        + for<'a> MulAssign<<Store::Scalar as Ref>::Ref<'a>>
        + 'static,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    FL: FunctionLibrary<Store::Tensor, Store::Scalar, Key = FK>,
    Aind: AbsInd,
    LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn execute_self_loop_traces(
        &mut self,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<bool, TensorNetworkError<K, FK>> {
        let mut did_trace = false;
        while execute_next_self_loop_trace::<LT, L, K, FK, Aind, Self>(self, graph, lib)? {
            did_trace = true;
            graph.sync_order();
        }

        Ok(did_trace)
    }

    fn execute_self_loop_traces_ignoring(
        &mut self,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        ignored: &mut SuBitGraph,
    ) -> Result<bool, TensorNetworkError<K, FK>> {
        let mut did_trace = false;
        while execute_next_self_loop_trace_ignoring::<LT, L, K, FK, Aind, Self>(
            self, graph, lib, ignored,
        )? {
            did_trace = true;
            graph.sync_order();
        }

        Ok(did_trace)
    }

    fn execute<C: ContractionStrategy<Self, L, K, FK, Aind>>(
        &mut self,
        graph: &NetworkGraph<K, FK, Aind>,
        operation: &NetworkOperation<FK>,
        lib: &L,
        fn_lib: &FL,
    ) -> Result<NetworkLeaf<K, Aind>, TensorNetworkError<K, FK>> {
        let _execute_span = profile::span(Timer::ExecuteOp);
        let op = operation.op().clone();
        let op_timer = match &op {
            NetworkOp::Neg => {
                profile::bump(Counter::ExecuteNeg, 1);
                Timer::ExecuteNeg
            }
            NetworkOp::Product => {
                profile::bump(Counter::ExecuteProduct, 1);
                Timer::ExecuteProduct
            }
            NetworkOp::Sum => {
                profile::bump(Counter::ExecuteSum, 1);
                Timer::ExecuteSum
            }
            NetworkOp::Function(_) => {
                profile::bump(Counter::ExecuteFunction, 1);
                Timer::ExecuteFunction
            }
            NetworkOp::Power(_) => {
                profile::bump(Counter::ExecutePower, 1);
                Timer::ExecutePower
            }
        };
        let _op_span = profile::span(op_timer);
        match op {
            NetworkOp::Neg => {
                if operation.children().len() != 1 {
                    return Err(TensorNetworkError::MoreThanOneNeg);
                }
                let child_id = operation.children()[0];
                if let NetworkNode::Leaf(leaf) = &graph.graph[child_id] {
                    let new_node = match leaf {
                        NetworkLeaf::Scalar(s) => {
                            let s = self.scalar_ref(*s).clone().neg();
                            let pos = self.push_scalar(s);

                            NetworkLeaf::Scalar(pos.into())
                        }
                        NetworkLeaf::LibraryKey { .. } => {
                            let inds = graph.get_lib_data(lib, child_id).unwrap();

                            let t = Store::Tensor::from(inds).neg();
                            let pos = self.push_tensor(t);
                            NetworkLeaf::LocalTensor(pos)
                        }
                        NetworkLeaf::LocalTensor(t) => {
                            let t = self.tensor(*t).clone().neg();
                            let pos = self.push_tensor(t);
                            NetworkLeaf::LocalTensor(pos)
                        }
                        NetworkLeaf::TensorSum(indices) => {
                            let mut negated = Vec::with_capacity(indices.len());
                            for index in indices {
                                let tensor = self.tensor(*index).clone().neg();
                                negated.push(self.push_tensor(tensor));
                            }
                            NetworkLeaf::TensorSum(negated)
                        }
                        NetworkLeaf::TensorTerm(term) => {
                            let term = if let Some(scalar) = term.scalar {
                                let scalar = self.scalar_ref(scalar).clone().neg();
                                TensorTerm::scaled(term.tensor, self.push_scalar(scalar))
                            } else {
                                let tensor = self.tensor(term.tensor).clone().neg();
                                TensorTerm::tensor(self.push_tensor(tensor))
                            };
                            NetworkLeaf::TensorTerm(term)
                        }
                        NetworkLeaf::TensorTermSum(terms) => {
                            let mut negated = Vec::with_capacity(terms.len());
                            for term in terms {
                                if let Some(scalar) = term.scalar {
                                    let scalar = self.scalar_ref(scalar).clone().neg();
                                    negated.push(TensorTerm::scaled(
                                        term.tensor,
                                        self.push_scalar(scalar),
                                    ));
                                } else {
                                    let tensor = self.tensor(term.tensor).clone().neg();
                                    negated.push(TensorTerm::tensor(self.push_tensor(tensor)));
                                }
                            }
                            NetworkLeaf::TensorTermSum(negated)
                        }
                    };
                    Ok(new_node)
                } else {
                    Err(TensorNetworkError::ChildlessNeg)
                }
            }
            NetworkOp::Product => {
                // println!("Doing Product");
                C::contract(self, graph, operation, lib)
            }
            NetworkOp::Sum => {
                let profile_sum = profile::enabled()
                    && (operation.leaf_count() > LARGE_SUM_PROFILE_THRESHOLD || profile::verbose());
                let sum_start = profile_sum.then(std::time::Instant::now);
                if profile_sum {
                    eprintln!(
                        "spenso_profile execute.sum_start leaves={} children={} subgraph_hedges={}",
                        operation.leaf_count(),
                        operation.children().len(),
                        operation.subgraph().n_included(),
                    );
                }

                let target_start = profile_sum.then(std::time::Instant::now);
                let mut targets = Vec::new();
                let mut scalar_targets = 0usize;
                let mut tensor_targets = 0usize;
                let mut library_targets = 0usize;
                for node in operation.children() {
                    if let NetworkNode::Leaf(l) = &graph.graph[*node] {
                        match l {
                            NetworkLeaf::Scalar(_) => scalar_targets += 1,
                            NetworkLeaf::LocalTensor(_)
                            | NetworkLeaf::TensorSum(_)
                            | NetworkLeaf::TensorTerm(_)
                            | NetworkLeaf::TensorTermSum(_) => tensor_targets += 1,
                            NetworkLeaf::LibraryKey { .. } => library_targets += 1,
                        }
                        targets.push((*node, l));
                    }
                }
                if let Some(target_start) = target_start {
                    eprintln!(
                        "spenso_profile execute.sum_targets leaves={} scalar_targets={} tensor_targets={} library_targets={} elapsed_ms={:.3}",
                        targets.len(),
                        scalar_targets,
                        tensor_targets,
                        library_targets,
                        target_start.elapsed().as_secs_f64() * 1000.0,
                    );
                }
                #[cfg(feature = "shadowing")]
                if profile_sum {
                    log_sum_leaf_atom_shapes::<K, FK, Aind, Self, LT, L>(
                        self, graph, &targets, lib,
                    );
                }

                let (nf, first) = &targets[0];

                #[cfg(feature = "shadowing")]
                if let Some(new_node) = try_atom_scalar_sum(self, &targets) {
                    if let Some(sum_start) = sum_start {
                        eprintln!(
                            "spenso_profile execute.sum_done leaves={} elapsed_ms={:.3}",
                            targets.len(),
                            sum_start.elapsed().as_secs_f64() * 1000.0,
                        );
                    }
                    return Ok(new_node);
                }

                if let Some(new_node) = try_balanced_scalar_sum(self, &targets, sum_start) {
                    if let Some(sum_start) = sum_start {
                        eprintln!(
                            "spenso_profile execute.sum_done leaves={} elapsed_ms={:.3}",
                            targets.len(),
                            sum_start.elapsed().as_secs_f64() * 1000.0,
                        );
                    }
                    return Ok(new_node);
                }

                if let Some(new_node) = try_balanced_tensor_sum::<K, FK, Aind, Self, LT, L>(
                    self, graph, &targets, lib, sum_start,
                ) {
                    if let Some(sum_start) = sum_start {
                        eprintln!(
                            "spenso_profile execute.sum_done leaves={} elapsed_ms={:.3}",
                            targets.len(),
                            sum_start.elapsed().as_secs_f64() * 1000.0,
                        );
                    }
                    return Ok(new_node);
                }

                let new_node = match first {
                    NetworkLeaf::Scalar(s) => {
                        let mut accumulator = self.scalar_ref(*s).clone();

                        for (offset, (_, t)) in targets[1..].iter().enumerate() {
                            let add_start = std::time::Instant::now();
                            match t {
                                NetworkLeaf::Scalar(s) => {
                                    accumulator += self.scalar_ref(*s).refer();
                                }
                                NetworkLeaf::LocalTensor(t) => {
                                    if let Some(s) = self.tensor(*t).scalar_ref() {
                                        accumulator += s;
                                    } else {
                                        return Err(TensorNetworkError::NotAllScalars(
                                            "".to_string(),
                                        ));
                                    }
                                }
                                NetworkLeaf::TensorSum(indices) => {
                                    for index in indices {
                                        if let Some(s) = self.tensor(*index).scalar_ref() {
                                            accumulator += s;
                                        } else {
                                            return Err(TensorNetworkError::NotAllScalars(
                                                "".to_string(),
                                            ));
                                        }
                                    }
                                }
                                NetworkLeaf::TensorTerm(term) => {
                                    let materialized = materialize_tensor_terms::<K, Aind, Self>(
                                        self,
                                        std::slice::from_ref(term),
                                        None,
                                    );
                                    let NetworkLeaf::LocalTensor(index) = materialized else {
                                        unreachable!("materialized tensor term is a local tensor")
                                    };
                                    if let Some(s) = self.tensor(index).scalar_ref() {
                                        accumulator += s;
                                    } else {
                                        return Err(TensorNetworkError::NotAllScalars(
                                            "".to_string(),
                                        ));
                                    }
                                }
                                NetworkLeaf::TensorTermSum(terms) => {
                                    let materialized = materialize_tensor_terms::<K, Aind, Self>(
                                        self, terms, None,
                                    );
                                    let NetworkLeaf::LocalTensor(index) = materialized else {
                                        unreachable!("materialized tensor terms are a local tensor")
                                    };
                                    if let Some(s) = self.tensor(index).scalar_ref() {
                                        accumulator += s;
                                    } else {
                                        return Err(TensorNetworkError::NotAllScalars(
                                            "".to_string(),
                                        ));
                                    }
                                }
                                NetworkLeaf::LibraryKey { .. } => {
                                    return Err(TensorNetworkError::ScalarLibSum("".to_string()));
                                }
                            }
                            log_sum_add(sum_start, offset + 2, targets.len(), add_start);
                        }

                        let pos = self.push_scalar(accumulator);
                        NetworkLeaf::Scalar(pos.into())
                    }
                    NetworkLeaf::LocalTensor(t) => {
                        let mut accumulator = self.tensor(*t).clone();
                        if accumulator.is_scalar() {
                            let mut accumulator =
                                Store::Scalar::from(accumulator.scalar().unwrap());

                            for (offset, (_, t)) in targets[1..].iter().enumerate() {
                                let add_start = std::time::Instant::now();
                                match t {
                                    NetworkLeaf::Scalar(s) => {
                                        accumulator += self.scalar_ref(*s).refer();
                                    }
                                    NetworkLeaf::LocalTensor(t) => {
                                        if let Some(s) = self.tensor(*t).scalar_ref() {
                                            accumulator += s;
                                        } else {
                                            return Err(TensorNetworkError::NotAllScalars(
                                                "".to_string(),
                                            ));
                                        }
                                    }
                                    NetworkLeaf::TensorSum(indices) => {
                                        for index in indices {
                                            if let Some(s) = self.tensor(*index).scalar_ref() {
                                                accumulator += s;
                                            } else {
                                                return Err(TensorNetworkError::NotAllScalars(
                                                    "".to_string(),
                                                ));
                                            }
                                        }
                                    }
                                    NetworkLeaf::TensorTerm(term) => {
                                        let materialized = materialize_tensor_terms::<K, Aind, Self>(
                                            self,
                                            std::slice::from_ref(term),
                                            None,
                                        );
                                        let NetworkLeaf::LocalTensor(index) = materialized else {
                                            unreachable!(
                                                "materialized tensor term is a local tensor"
                                            )
                                        };
                                        if let Some(s) = self.tensor(index).scalar_ref() {
                                            accumulator += s;
                                        } else {
                                            return Err(TensorNetworkError::NotAllScalars(
                                                "".to_string(),
                                            ));
                                        }
                                    }
                                    NetworkLeaf::TensorTermSum(terms) => {
                                        let materialized = materialize_tensor_terms::<K, Aind, Self>(
                                            self, terms, None,
                                        );
                                        let NetworkLeaf::LocalTensor(index) = materialized else {
                                            unreachable!(
                                                "materialized tensor terms are a local tensor"
                                            )
                                        };
                                        if let Some(s) = self.tensor(index).scalar_ref() {
                                            accumulator += s;
                                        } else {
                                            return Err(TensorNetworkError::NotAllScalars(
                                                "".to_string(),
                                            ));
                                        }
                                    }
                                    NetworkLeaf::LibraryKey { .. } => {
                                        return Err(TensorNetworkError::ScalarLibSum(
                                            "".to_string(),
                                        ));
                                    }
                                }
                                log_sum_add(sum_start, offset + 2, targets.len(), add_start);
                            }

                            let pos = self.push_scalar(accumulator);
                            NetworkLeaf::Scalar(pos.into())
                        } else {
                            for (offset, (nid, t)) in targets[1..].iter().enumerate() {
                                let add_start = std::time::Instant::now();
                                match t {
                                    NetworkLeaf::Scalar(_) => {
                                        return Err(TensorNetworkError::SumScalarTensor(
                                            "".to_string(),
                                        ));
                                    }
                                    NetworkLeaf::LocalTensor(t) => {
                                        accumulator += self.tensor(*t).refer();
                                    }
                                    NetworkLeaf::TensorSum(indices) => {
                                        for index in indices {
                                            accumulator += self.tensor(*index).refer();
                                        }
                                    }
                                    NetworkLeaf::TensorTerm(term) => {
                                        let materialized = materialize_tensor_terms::<K, Aind, Self>(
                                            self,
                                            std::slice::from_ref(term),
                                            None,
                                        );
                                        let NetworkLeaf::LocalTensor(index) = materialized else {
                                            unreachable!(
                                                "materialized tensor term is a local tensor"
                                            )
                                        };
                                        accumulator += self.tensor(index).refer();
                                    }
                                    NetworkLeaf::TensorTermSum(terms) => {
                                        let materialized = materialize_tensor_terms::<K, Aind, Self>(
                                            self, terms, None,
                                        );
                                        let NetworkLeaf::LocalTensor(index) = materialized else {
                                            unreachable!(
                                                "materialized tensor terms are a local tensor"
                                            )
                                        };
                                        accumulator += self.tensor(index).refer();
                                    }
                                    NetworkLeaf::LibraryKey { .. } => {
                                        let with_index = graph.get_lib_data(lib, *nid).unwrap();

                                        accumulator += with_index;
                                    }
                                }
                                log_sum_add(sum_start, offset + 2, targets.len(), add_start);
                            }

                            let pos = self.push_tensor(accumulator);

                            NetworkLeaf::LocalTensor(pos)
                        }
                    }
                    NetworkLeaf::LibraryKey { .. } => {
                        let inds = graph.get_lib_data(lib, *nf).unwrap();
                        let mut accumulator = Store::Tensor::from(inds);
                        for (offset, (nid, t)) in targets[1..].iter().enumerate() {
                            let add_start = std::time::Instant::now();
                            match t {
                                NetworkLeaf::Scalar(_) => {
                                    return Err(TensorNetworkError::SumScalarTensor(
                                        "".to_string(),
                                    ));
                                }
                                NetworkLeaf::LocalTensor(t) => {
                                    accumulator += self.tensor(*t).refer();
                                }
                                NetworkLeaf::TensorSum(indices) => {
                                    for index in indices {
                                        accumulator += self.tensor(*index).refer();
                                    }
                                }
                                NetworkLeaf::TensorTerm(term) => {
                                    let materialized = materialize_tensor_terms::<K, Aind, Self>(
                                        self,
                                        std::slice::from_ref(term),
                                        None,
                                    );
                                    let NetworkLeaf::LocalTensor(index) = materialized else {
                                        unreachable!("materialized tensor term is a local tensor")
                                    };
                                    accumulator += self.tensor(index).refer();
                                }
                                NetworkLeaf::TensorTermSum(terms) => {
                                    let materialized = materialize_tensor_terms::<K, Aind, Self>(
                                        self, terms, None,
                                    );
                                    let NetworkLeaf::LocalTensor(index) = materialized else {
                                        unreachable!("materialized tensor terms are a local tensor")
                                    };
                                    accumulator += self.tensor(index).refer();
                                }
                                NetworkLeaf::LibraryKey { .. } => {
                                    let with = graph.get_lib_data(lib, *nid).unwrap();
                                    accumulator += with;
                                }
                            }
                            log_sum_add(sum_start, offset + 2, targets.len(), add_start);
                        }

                        let pos = self.push_tensor(accumulator);

                        NetworkLeaf::LocalTensor(pos)
                    }
                    NetworkLeaf::TensorSum(_)
                    | NetworkLeaf::TensorTerm(_)
                    | NetworkLeaf::TensorTermSum(_) => {
                        return Err(TensorNetworkError::SumScalarTensor(
                            "lazy tensor sum mixed with scalar sum terms".to_string(),
                        ));
                    }
                };
                if let Some(sum_start) = sum_start {
                    eprintln!(
                        "spenso_profile execute.sum_done leaves={} elapsed_ms={:.3}",
                        targets.len(),
                        sum_start.elapsed().as_secs_f64() * 1000.0,
                    );
                }

                Ok(new_node)
            }
            NetworkOp::Function(f) => {
                if operation.children().len() != 1 {
                    return Err(TensorNetworkError::Other(eyre!(
                        "Cannot have more than one tensor argument to function"
                    )));
                }
                let child_id = operation.children()[0];
                if let NetworkNode::Leaf(leaf) = &graph.graph[child_id] {
                    let new_node = match leaf {
                        NetworkLeaf::Scalar(s) => {
                            let s = self.scalar_ref(*s).clone();
                            let s = fn_lib.apply_scalar(&f, s)?;
                            let pos = self.push_scalar(s);

                            NetworkLeaf::Scalar(pos.into())
                        }
                        NetworkLeaf::LibraryKey { .. } => {
                            let inds = graph.get_lib_data(lib, child_id).unwrap();
                            let t = fn_lib.apply(&f, Store::Tensor::from(inds))?;
                            let pos = self.push_tensor(t);
                            NetworkLeaf::LocalTensor(pos)
                        }
                        NetworkLeaf::LocalTensor(t) => {
                            let t = self.tensor(*t).clone();
                            let t = fn_lib.apply(&f, t)?;
                            let pos = self.push_tensor(t);
                            NetworkLeaf::LocalTensor(pos)
                        }
                        NetworkLeaf::TensorSum(indices) => {
                            let materialized =
                                materialize_tensor_sum::<K, Aind, Self>(self, indices, None);
                            let NetworkLeaf::LocalTensor(index) = materialized else {
                                unreachable!("materialized tensor sum is a local tensor")
                            };
                            let t = self.tensor(index).clone();
                            let t = fn_lib.apply(&f, t)?;
                            let pos = self.push_tensor(t);
                            NetworkLeaf::LocalTensor(pos)
                        }
                        NetworkLeaf::TensorTerm(term) => {
                            let materialized = materialize_tensor_terms::<K, Aind, Self>(
                                self,
                                std::slice::from_ref(term),
                                None,
                            );
                            let NetworkLeaf::LocalTensor(index) = materialized else {
                                unreachable!("materialized tensor term is a local tensor")
                            };
                            let t = self.tensor(index).clone();
                            let t = fn_lib.apply(&f, t)?;
                            let pos = self.push_tensor(t);
                            NetworkLeaf::LocalTensor(pos)
                        }
                        NetworkLeaf::TensorTermSum(terms) => {
                            let materialized =
                                materialize_tensor_terms::<K, Aind, Self>(self, terms, None);
                            let NetworkLeaf::LocalTensor(index) = materialized else {
                                unreachable!("materialized tensor terms are a local tensor")
                            };
                            let t = self.tensor(index).clone();
                            let t = fn_lib.apply(&f, t)?;
                            let pos = self.push_tensor(t);
                            NetworkLeaf::LocalTensor(pos)
                        }
                    };
                    Ok(new_node)
                } else {
                    Err(TensorNetworkError::ChildlessNeg)
                }
            }
            NetworkOp::Power(pow) => {
                if operation.children().len() != 1 {
                    return Err(TensorNetworkError::Other(eyre!(
                        "Cannot have more than one tensor argument to power:{}",
                        graph.dot()
                    )));
                }
                let n = pow.abs();
                let child_id = operation.children()[0];
                if let NetworkNode::Leaf(leaf) = &graph.graph[child_id] {
                    let new_node = match leaf {
                        NetworkLeaf::Scalar(si) => {
                            if n == 0 {
                                NetworkLeaf::Scalar(*si)
                            } else {
                                let mut s = self.scalar_ref(*si).clone();

                                for _ in 1..n {
                                    s *= self.scalar_ref(*si).refer();
                                }

                                if pow < 0 {
                                    s = s.ref_one() / s;
                                }

                                let pos = self.push_scalar(s);

                                NetworkLeaf::Scalar(pos.into())
                            }
                        }
                        NetworkLeaf::LibraryKey { key, indices } => {
                            let inds = graph.get_lib_data(lib, child_id).unwrap();
                            let mut t = Store::Tensor::from(inds);

                            match pow {
                                0 => {
                                    let one = self.scalar(0).ref_one();
                                    let pos = self.push_scalar(one);
                                    NetworkLeaf::Scalar(pos.into())
                                }
                                1 => NetworkLeaf::LibraryKey {
                                    key: key.clone(),
                                    indices: indices.clone(),
                                },
                                _ => {
                                    let squares = n / 2;
                                    let mut square = t.contract(&t)?;

                                    if n % 2 == 1 {
                                        if n != 1 {
                                            for _ in 0..squares {
                                                square = square.contract(&square)?;
                                            }
                                            t = square.contract(&t)?;
                                        }

                                        if pow < 0 {
                                            if !t.is_scalar() {
                                                return Err(
                                                    TensorNetworkError::NegativeExponentNonScalar(
                                                        "".to_string(),
                                                    ),
                                                );
                                            } else {
                                                let mut s =
                                                    Store::Scalar::from(t.scalar().unwrap());
                                                s = s.ref_one() / s;
                                                let pos = self.push_scalar(s);
                                                NetworkLeaf::Scalar(pos.into())
                                            }
                                        } else {
                                            let pos = self.push_tensor(t);
                                            NetworkLeaf::LocalTensor(pos)
                                        }
                                    } else {
                                        let mut s = Store::Scalar::from(square.scalar().unwrap());
                                        let sc = s.clone();
                                        for _ in 1..squares {
                                            s *= sc.refer();
                                        }
                                        if pow < 0 {
                                            s = s.ref_one() / s;
                                        }
                                        let pos = self.push_scalar(s);
                                        NetworkLeaf::Scalar(pos.into())
                                    }
                                }
                            }
                        }
                        NetworkLeaf::LocalTensor(ti) => {
                            let mut t = self.tensor(*ti).clone();
                            match pow {
                                0 => {
                                    let one = self.scalar(0).ref_one();
                                    let pos = self.push_scalar(one);
                                    NetworkLeaf::Scalar(pos.into())
                                }
                                1 => NetworkLeaf::LocalTensor(*ti),
                                _ => {
                                    let squares = n / 2;
                                    let mut square = t.contract(&t)?;

                                    if n % 2 == 1 {
                                        if n != 1 {
                                            for _ in 0..squares {
                                                square = square.contract(&square)?;
                                            }
                                            t = square.contract(&t)?;
                                        }
                                        if pow < 0 {
                                            if !t.is_scalar() {
                                                return Err(
                                                    TensorNetworkError::NegativeExponentNonScalar(
                                                        "".to_string(),
                                                    ),
                                                );
                                            } else {
                                                let mut s =
                                                    Store::Scalar::from(t.scalar().unwrap());
                                                s = s.ref_one() / s;
                                                let pos = self.push_scalar(s);
                                                NetworkLeaf::Scalar(pos.into())
                                            }
                                        } else {
                                            let pos = self.push_tensor(t);
                                            NetworkLeaf::LocalTensor(pos)
                                        }
                                    } else {
                                        let mut s = Store::Scalar::from(square.scalar().unwrap());
                                        let sc = s.clone();
                                        for _ in 1..squares {
                                            s *= sc.refer();
                                        }
                                        if pow < 0 {
                                            s = s.ref_one() / s;
                                        }
                                        let pos = self.push_scalar(s);
                                        NetworkLeaf::Scalar(pos.into())
                                    }
                                }
                            }
                        }
                        NetworkLeaf::TensorSum(indices) => {
                            let materialized =
                                materialize_tensor_sum::<K, Aind, Self>(self, indices, None);
                            let NetworkLeaf::LocalTensor(ti) = materialized else {
                                unreachable!("materialized tensor sum is a local tensor")
                            };
                            let mut t = self.tensor(ti).clone();
                            match pow {
                                0 => {
                                    let one = self.scalar(0).ref_one();
                                    let pos = self.push_scalar(one);
                                    NetworkLeaf::Scalar(pos.into())
                                }
                                1 => NetworkLeaf::LocalTensor(ti),
                                _ => {
                                    let squares = n / 2;
                                    let mut square = t.contract(&t)?;

                                    if n % 2 == 1 {
                                        if n != 1 {
                                            for _ in 0..squares {
                                                square = square.contract(&square)?;
                                            }
                                            t = square.contract(&t)?;
                                        }
                                        if pow < 0 {
                                            if !t.is_scalar() {
                                                return Err(
                                                    TensorNetworkError::NegativeExponentNonScalar(
                                                        "".to_string(),
                                                    ),
                                                );
                                            } else {
                                                let mut s =
                                                    Store::Scalar::from(t.scalar().unwrap());
                                                s = s.ref_one() / s;
                                                let pos = self.push_scalar(s);
                                                NetworkLeaf::Scalar(pos.into())
                                            }
                                        } else {
                                            let pos = self.push_tensor(t);
                                            NetworkLeaf::LocalTensor(pos)
                                        }
                                    } else {
                                        let mut s = Store::Scalar::from(square.scalar().unwrap());
                                        let sc = s.clone();
                                        for _ in 1..squares {
                                            s *= sc.refer();
                                        }
                                        if pow < 0 {
                                            s = s.ref_one() / s;
                                        }
                                        let pos = self.push_scalar(s);
                                        NetworkLeaf::Scalar(pos.into())
                                    }
                                }
                            }
                        }
                        NetworkLeaf::TensorTerm(term) => {
                            let materialized = materialize_tensor_terms::<K, Aind, Self>(
                                self,
                                std::slice::from_ref(term),
                                None,
                            );
                            let NetworkLeaf::LocalTensor(ti) = materialized else {
                                unreachable!("materialized tensor term is a local tensor")
                            };
                            let mut t = self.tensor(ti).clone();
                            match pow {
                                0 => {
                                    let one = self.scalar(0).ref_one();
                                    let pos = self.push_scalar(one);
                                    NetworkLeaf::Scalar(pos.into())
                                }
                                1 => NetworkLeaf::LocalTensor(ti),
                                _ => {
                                    let squares = n / 2;
                                    let mut square = t.contract(&t)?;

                                    if n % 2 == 1 {
                                        if n != 1 {
                                            for _ in 0..squares {
                                                square = square.contract(&square)?;
                                            }
                                            t = square.contract(&t)?;
                                        }
                                        if pow < 0 {
                                            if !t.is_scalar() {
                                                return Err(
                                                    TensorNetworkError::NegativeExponentNonScalar(
                                                        "".to_string(),
                                                    ),
                                                );
                                            } else {
                                                let mut s =
                                                    Store::Scalar::from(t.scalar().unwrap());
                                                s = s.ref_one() / s;
                                                let pos = self.push_scalar(s);
                                                NetworkLeaf::Scalar(pos.into())
                                            }
                                        } else {
                                            let pos = self.push_tensor(t);
                                            NetworkLeaf::LocalTensor(pos)
                                        }
                                    } else {
                                        let mut s = Store::Scalar::from(square.scalar().unwrap());
                                        let sc = s.clone();
                                        for _ in 1..squares {
                                            s *= sc.refer();
                                        }
                                        if pow < 0 {
                                            s = s.ref_one() / s;
                                        }
                                        let pos = self.push_scalar(s);
                                        NetworkLeaf::Scalar(pos.into())
                                    }
                                }
                            }
                        }
                        NetworkLeaf::TensorTermSum(terms) => {
                            let materialized =
                                materialize_tensor_terms::<K, Aind, Self>(self, terms, None);
                            let NetworkLeaf::LocalTensor(ti) = materialized else {
                                unreachable!("materialized tensor terms are a local tensor")
                            };
                            let mut t = self.tensor(ti).clone();
                            match pow {
                                0 => {
                                    let one = self.scalar(0).ref_one();
                                    let pos = self.push_scalar(one);
                                    NetworkLeaf::Scalar(pos.into())
                                }
                                1 => NetworkLeaf::LocalTensor(ti),
                                _ => {
                                    let squares = n / 2;
                                    let mut square = t.contract(&t)?;

                                    if n % 2 == 1 {
                                        if n != 1 {
                                            for _ in 0..squares {
                                                square = square.contract(&square)?;
                                            }
                                            t = square.contract(&t)?;
                                        }
                                        if pow < 0 {
                                            if !t.is_scalar() {
                                                return Err(
                                                    TensorNetworkError::NegativeExponentNonScalar(
                                                        "".to_string(),
                                                    ),
                                                );
                                            } else {
                                                let mut s =
                                                    Store::Scalar::from(t.scalar().unwrap());
                                                s = s.ref_one() / s;
                                                let pos = self.push_scalar(s);
                                                NetworkLeaf::Scalar(pos.into())
                                            }
                                        } else {
                                            let pos = self.push_tensor(t);
                                            NetworkLeaf::LocalTensor(pos)
                                        }
                                    } else {
                                        let mut s = Store::Scalar::from(square.scalar().unwrap());
                                        let sc = s.clone();
                                        for _ in 1..squares {
                                            s *= sc.refer();
                                        }
                                        if pow < 0 {
                                            s = s.ref_one() / s;
                                        }
                                        let pos = self.push_scalar(s);
                                        NetworkLeaf::Scalar(pos.into())
                                    }
                                }
                            }
                        }
                    };
                    Ok(new_node)
                } else {
                    Err(TensorNetworkError::ChildlessNeg)
                }
            }
        }
    }
}

pub trait Ref {
    type Ref<'a>
    where
        Self: 'a;
    fn refer(&self) -> Self::Ref<'_>;
}

impl Ref for f64 {
    type Ref<'a>
        = &'a f64
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        self
    }
}

// #[cfg(feature = "shadowing")]
// pub mod levels;
#[cfg(feature = "shadowing")]
pub mod symbolica_interop;

#[cfg(test)]
mod tests;
