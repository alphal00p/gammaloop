use graph::{
    NAdd, NMul, NetworkEdge, NetworkGraph, NetworkLeaf, NetworkNode, NetworkOp, NetworkOperation,
};
use linnet::half_edge::{
    NodeIndex,
    subgraph::{SuBitGraph, SubSetLike},
};
use profile::{Counter, Timer};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use library::{Library, LibraryError};

use crate::algebra::algebraic_traits::RefOne;
use crate::contraction::{Contract, Trace};
use crate::network::library::{DummyKey, FunctionLibrary, FunctionLibraryError, LibraryTensor};
use crate::structure::abstract_index::AbstractIndex;
use crate::structure::permuted::PermuteTensor;
// use crate::shadowing::Concretize;
use crate::structure::representation::LibrarySlot;
use crate::structure::slot::{AbsInd, IsAbstractSlot};
use crate::structure::{HasName, PermutedStructure, StructureError};
use std::borrow::Cow;
use std::fmt::Display;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};
use store::{
    NetworkStore, NetworkStoreAccess, NetworkStoreOverlay, TensorScalarStore,
    TensorScalarStoreMapping,
};
use thiserror::Error;
// use log::trace;

use crate::{
    contraction::ContractionError,
    structure::{CastStructure, HasStructure, ScalarTensor, TensorStructure},
};
use eyre::eyre;

#[cfg(feature = "shadowing")]
pub mod tags;
// use eyre::Result;

use std::{convert::Infallible, fmt::Debug};

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
                NetworkLeaf::Scalar(t) => Ok(ExecutionResult::Val(TensorOrScalarOrKey::Scalar(
                    self.store.get_scalar(*t),
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
        T: Clone + ScalarTensor + HasStructure,
        K: Display + Debug,
        FK: Display + Debug + Clone,
        LT: TensorStructure<Indexed = T> + Clone + LibraryTensor<WithIndices = T>,
        T: PermuteTensor<Permuted = T>,
        for<'b> &'b S: Into<T::Scalar>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        Ok(match self.result()? {
            ExecutionResult::One => ExecutionResult::One,
            ExecutionResult::Zero => ExecutionResult::Zero,
            ExecutionResult::Val(v) => ExecutionResult::Val(match v {
                TensorOrScalarOrKey::Tensor { tensor, .. } => Cow::Borrowed(tensor),
                TensorOrScalarOrKey::Scalar(s) => Cow::Owned(T::new_scalar(s.into())),
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
        S: Clone,
    {
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
                    NetworkLeaf::Scalar(s) => Some(format!(
                        "label = \"S:{}\"",
                        scalar_disp(self.store.get_scalar(*s))
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
    ContractScalars, ContractionStrategy, MinResultRank, ProductContraction, SingleSmallestDegree,
    SmallestDegree, SmallestDegreeIter,
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

pub struct Sequential;
pub struct Parallel;
pub struct SequentialRef;
pub struct SequentialExtract;

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
        let mut batch_index = 0usize;

        loop {
            if executor.execute_self_loop_traces(graph, lib)? {
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
}

fn remap_parallel_replacement<K, Aind>(
    replacement: &mut NetworkLeaf<K, Aind>,
    base_tensors: usize,
    base_scalars: usize,
    tensor_offset: usize,
    scalar_offset: usize,
) {
    match replacement {
        NetworkLeaf::LocalTensor(index) if *index >= base_tensors => {
            *index = tensor_offset + (*index - base_tensors);
        }
        NetworkLeaf::Scalar(index) if *index >= base_scalars => {
            *index = scalar_offset + (*index - base_scalars);
        }
        NetworkLeaf::LocalTensor(_) | NetworkLeaf::LibraryKey { .. } | NetworkLeaf::Scalar(_) => {}
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
                        let (tensors, scalars) = local.into_additions();

                        Ok(ParallelExecutionOutput {
                            replacement,
                            tensors,
                            scalars,
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

impl<T, Sc> NetworkStore<T, Sc> {
    fn execute_next_self_loop_trace<LT, L, K, FK, Aind>(
        &mut self,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<bool, TensorNetworkError<K, FK>>
    where
        LT: LibraryTensor + Clone,
        T: HasStructure + TensorStructure + Clone + Trace + From<LT::WithIndices>,
        L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
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
                self.tensors[*local].internal_contract()
            }
            NetworkNode::Leaf(NetworkLeaf::LibraryKey { .. }) => T::from(
                graph
                    .get_lib_data(lib, node)
                    .expect("library node selected for trace must have library data"),
            )
            .internal_contract(),
            NetworkNode::Leaf(NetworkLeaf::Scalar(_)) => return Ok(false),
            NetworkNode::Op(_) => unreachable!("self-loop trace search only selects tensor leaves"),
        };

        let traced_position = self.tensors.len();
        self.tensors.push(traced);
        graph.replace_node_deleting_self_loop_slots(
            node,
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(traced_position)),
        );

        Ok(true)
    }
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
        + Contract<LCM = Store::Tensor>
        + for<'a> AddAssign<<Store::Tensor as Ref>::Ref<'a>>
        + for<'a> AddAssign<LT::WithIndices>
        + From<LT::WithIndices>,
    L: Library<<Store::Tensor as HasStructure>::Structure, Key = K, Value = PermutedStructure<LT>>,
    Store::Scalar: Neg<Output = Store::Scalar>
        + RefOne
        + Div<Output = Store::Scalar>
        + for<'a> AddAssign<<Store::Scalar as Ref>::Ref<'a>>
        + Clone
        + for<'a> AddAssign<<Store::Tensor as HasStructure>::ScalarRef<'a>>
        + From<<Store::Tensor as HasStructure>::Scalar>
        + Ref
        + for<'a> MulAssign<<Store::Scalar as Ref>::Ref<'a>>,
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
        while self.execute_next_self_loop_trace(graph, lib)? {
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
                            let s = self.scalar(*s).clone().neg();
                            let pos = self.push_scalar(s);

                            NetworkLeaf::Scalar(pos)
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
                let profile_sum = profile::enabled() && operation.leaf_count() > 32;
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
                            NetworkLeaf::LocalTensor(_) => tensor_targets += 1,
                            NetworkLeaf::LibraryKey(_) => library_targets += 1,
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

                let log_sum_add = |done: usize, total: usize, add_start: std::time::Instant| {
                    if let Some(sum_start) = sum_start {
                        let add_elapsed = add_start.elapsed();
                        if add_elapsed.as_millis() >= 1_000 || done <= 8 || done.is_multiple_of(16)
                        {
                            eprintln!(
                                "spenso_profile execute.sum_progress done={} total={} add_ms={:.3} elapsed_ms={:.3}",
                                done,
                                total,
                                add_elapsed.as_secs_f64() * 1000.0,
                                sum_start.elapsed().as_secs_f64() * 1000.0,
                            );
                        }
                    }
                };

                let (nf, first) = &targets[0];

                let new_node = match first {
                    NetworkLeaf::Scalar(s) => {
                        let mut accumulator = self.scalar(*s).clone();

                        for (offset, (_, t)) in targets[1..].iter().enumerate() {
                            let add_start = std::time::Instant::now();
                            match t {
                                NetworkLeaf::Scalar(s) => {
                                    accumulator += self.scalar(*s).refer();
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
                                NetworkLeaf::LibraryKey { .. } => {
                                    return Err(TensorNetworkError::ScalarLibSum("".to_string()));
                                }
                            }
                            log_sum_add(offset + 2, targets.len(), add_start);
                        }

                        let pos = self.push_scalar(accumulator);
                        NetworkLeaf::Scalar(pos)
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
                                        accumulator += self.scalar(*s).refer();
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
                                    NetworkLeaf::LibraryKey { .. } => {
                                        return Err(TensorNetworkError::ScalarLibSum(
                                            "".to_string(),
                                        ));
                                    }
                                }
                                log_sum_add(offset + 2, targets.len(), add_start);
                            }

                            let pos = self.push_scalar(accumulator);
                            NetworkLeaf::Scalar(pos)
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
                                    NetworkLeaf::LibraryKey { .. } => {
                                        let with_index = graph.get_lib_data(lib, *nid).unwrap();

                                        accumulator += with_index;
                                    }
                                }
                                log_sum_add(offset + 2, targets.len(), add_start);
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
                                NetworkLeaf::LibraryKey { .. } => {
                                    let with = graph.get_lib_data(lib, *nid).unwrap();
                                    accumulator += with;
                                }
                            }
                            log_sum_add(offset + 2, targets.len(), add_start);
                        }

                        let pos = self.push_tensor(accumulator);

                        NetworkLeaf::LocalTensor(pos)
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
                            let s = self.scalar(*s).clone();
                            let s = fn_lib.apply_scalar(&f, s)?;
                            let pos = self.push_scalar(s);

                            NetworkLeaf::Scalar(pos)
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
                                let mut s = self.scalar(*si).clone();

                                for _ in 1..n {
                                    s *= self.scalar(*si).refer();
                                }

                                if pow < 0 {
                                    s = s.ref_one() / s;
                                }

                                let pos = self.push_scalar(s);

                                NetworkLeaf::Scalar(pos)
                            }
                        }
                        NetworkLeaf::LibraryKey { key, indices } => {
                            let inds = graph.get_lib_data(lib, child_id).unwrap();
                            let mut t = Store::Tensor::from(inds);

                            match pow {
                                0 => {
                                    let one = self.scalar(0).ref_one();
                                    let pos = self.push_scalar(one);
                                    NetworkLeaf::Scalar(pos)
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
                                                NetworkLeaf::Scalar(pos)
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
                                        NetworkLeaf::Scalar(pos)
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
                                    NetworkLeaf::Scalar(pos)
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
                                                NetworkLeaf::Scalar(pos)
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
                                        NetworkLeaf::Scalar(pos)
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
